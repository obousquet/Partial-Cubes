/* Horizontal-only refinement and progressive disclosure for ranked graphs.
 * No inferred edges or ranks are introduced here. The crossing score uses
 * chords of the displayed edges, a proxy for the renderer's gently bent curves.
 */
window.GraphReadability = (() => {
    const cross = (a, b, c) => (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
    function metrics(nodes, edges) {
        const byId = new Map(nodes.map(n => [n.id, n]));
        const lines = edges.map(e => [byId.get(e.source), byId.get(e.target)])
            .filter(([a,b]) => a && b && a !== b);
        let crossings = 0, span = 0;
        for (let i = 0; i < lines.length; i++) {
            const [a,b] = lines[i];
            span += Math.abs(a.x-b.x);
            for (let j = 0; j < i; j++) {
                const [c,d] = lines[j];
                if (a===c || a===d || b===c || b===d) continue;
                if (cross(a,b,c)*cross(a,b,d) < 0 && cross(c,d,a)*cross(c,d,b) < 0) crossings++;
            }
        }
        const width = Math.max(...nodes.map(n=>n.x+n.width/2)) - Math.min(...nodes.map(n=>n.x-n.width/2));
        return {crossings, span, width, score: crossings * 1000 + span / 20 + width / 5};
    }
    function optimize(nodes, edges, clusters) {
        if (!nodes.length) return {before: null, after: null};
        const unconstrained = metrics(nodes, edges), byId = new Map(nodes.map(n=>[n.id,n]));
        const grouped=nodes.every(n=>Number.isInteger(n.horizontal_group));
        const rows = new Map(), membership = new Map();
        clusters.forEach((c,i) => c.nodes.forEach(id=>membership.set(id, 'cluster-'+i)));
        for (const n of nodes) {
            const key = `${n.rank ?? n.y}:${grouped?n.horizontal_group:''}`;
            if (!rows.has(key)) rows.set(key, new Map());
            const units = rows.get(key), unit = membership.get(n.id) ?? n.id;
            if (!units.has(unit)) units.set(unit, []);
            units.get(unit).push(n);
        }
        const layers = [...rows.values()].map(units => [...units.values()]
            .map(unit=>unit.sort((a,b)=>a.x-b.x)).sort((a,b)=>a[0].x-b[0].x))
            .sort((a,b)=>a[0][0].y-b[0][0].y);
        const unitWidth=unit=>unit.reduce((s,n)=>s+n.width,0)+24*(unit.length-1)+32;
        const rowWidth=row=>row.reduce((s,unit)=>s+unitWidth(unit),0)+32*(row.length-1);
        const bands=new Map();
        if (grouped) {
            // Regions are global across rows, not merely a sort within each
            // rank. Allocate only the width required by the widest row of a
            // region; no dummy nodes or invisible routing obstacles.
            for(const row of layers) {
                const group=row[0][0].horizontal_group;
                bands.set(group,Math.max(bands.get(group)||0,rowWidth(row)));
            }
            let left=0;
            for(const group of [...bands.keys()].sort((a,b)=>a-b)) {
                const width=bands.get(group);
                bands.set(group,{left,right:left+width});
                left+=width+100;
            }
            for(const row of layers) {
                const band=bands.get(row[0][0].horizontal_group);
                let x=(band.left+band.right-rowWidth(row))/2+16;
                for(const unit of row) {
                    for(const n of unit) {n.x=x+n.width/2;x+=n.width+24;}
                    x+=40;
                }
            }
        }
        const withinBand=members=>!grouped||members.every(n=> {
            const band=bands.get(n.horizontal_group);
            return n.x-n.width/2>=band.left && n.x+n.width/2<=band.right;
        });
        const before=metrics(nodes,edges);
        const adjacency = new Map(nodes.map(n=>[n.id,[]]));
        edges.forEach(e=> {
            if (byId.has(e.source) && byId.has(e.target)) {
                adjacency.get(e.source).push(byId.get(e.target));
                adjacency.get(e.target).push(byId.get(e.source));
            }
        });
        let best = before.score;
        function tryRow(order, compact=true) {
            const members = order.flat(), saved = members.map(n=>n.x);
            const outside = members.flatMap(n=>adjacency.get(n.id)).filter(n=>!members.includes(n));
            let center = outside.length ? outside.reduce((s,n)=>s+n.x,0)/outside.length
                : saved.reduce((a,b)=>a+b,0)/saved.length;
            // Keep affine-linear blocks together with room for their frame.
            const widths = order.map(unitWidth);
            if(grouped) {
                const band=bands.get(members[0].horizontal_group), half=rowWidth(order)/2;
                center=Math.max(band.left+half,Math.min(band.right-half,center));
            }
            let left = center - (widths.reduce((a,b)=>a+b,0)+32*(order.length-1))/2;
            order.forEach((unit,i)=> {
                let x = left+16;
                unit.forEach(n=> { n.x=x+n.width/2; x+=n.width+24; });
                left+=widths[i]+32;
            });
            if (!compact) {
                const slots=[...saved].sort((a,b)=>a-b);
                let cursor=-Infinity, slot=0;
                order.forEach(unit=> {
                    unit.forEach(n=> {
                        n.x=Math.max(slots[slot++],cursor+n.width/2+32);
                        cursor=n.x+n.width/2;
                    });
                });
            }
            const candidate = metrics(nodes,edges).score;
            if (withinBand(members) && candidate < best - 1e-6) { best=candidate; return true; }
            members.forEach((n,i)=>n.x=saved[i]);
            return false;
        }
        // Test compact and existing-spacing candidates, rather than forcing
        // compactness at the expense of hundreds of new crossings.
        // Alternating barycentre sweeps, then bounded adjacent swaps. Only
        // accept global improvements, so overlays count as well as backbone.
        for (let pass=0; pass<6; pass++) {
            const sweep = pass%2 ? [...layers].reverse() : layers;
            for (const row of sweep) {
                const barycentre = unit => {
                    const neighbours=unit.flatMap(n=>adjacency.get(n.id));
                    return (neighbours.length?neighbours:unit).reduce((s,n)=>s+n.x,0)/(neighbours.length||unit.length);
                };
                const ordered=[...row].sort((a,b)=>barycentre(a)-barycentre(b));
                if (tryRow(ordered,false)) row.splice(0,row.length,...ordered);
                if (tryRow(ordered)) row.splice(0,row.length,...ordered);
                tryRow(row);
                for (let i=0;i+1<row.length;i++) {
                    const swapped=[...row]; [swapped[i],swapped[i+1]]=[swapped[i+1],swapped[i]];
                    if (tryRow(swapped,false)) row.splice(0,row.length,...swapped);
                    if (tryRow(swapped)) row.splice(0,row.length,...swapped);
                }
                // Move blocks toward their neighbours without changing the
                // row order. Preserve a hard clearance from adjacent blocks.
                for (let i=0;i<row.length;i++) {
                    const unit=row[i], saved=unit.map(n=>n.x);
                    const left=unit[0].x-unit[0].width/2;
                    const last=unit[unit.length-1], right=last.x+last.width/2;
                    const previous=i?row[i-1].at(-1):null, next=row[i+1]?.[0];
                    const band=grouped?bands.get(unit[0].horizontal_group):null;
                    const lo=previous?previous.x+previous.width/2+48-left:band?band.left+16-left:-before.width;
                    const hi=next?next.x-next.width/2-48-right:band?band.right-16-right:before.width;
                    if (lo>hi) continue;
                    const desired=barycentre(unit)-(left+right)/2;
                    for(const fraction of [1,0.5,0.25]) {
                        const shift=Math.max(lo,Math.min(hi,desired*fraction));
                        unit.forEach((n,j)=>n.x=saved[j]+shift);
                        const candidate=metrics(nodes,edges).score;
                        if(candidate<best-1e-6) {best=candidate;break;}
                        unit.forEach((n,j)=>n.x=saved[j]);
                    }
                }
            }
        }
        return {unconstrained, before, after: metrics(nodes,edges), grouped};
    }
    function place(svg, metadata, edges, clusters, enabled) {
        if (!svg || !enabled) return;
        const meta = new Map(metadata.map(n=>[n.id,n]));
        const nodes = [...svg.querySelectorAll('g.node')].map(element=> {
            const id=element.querySelector('title').textContent, box=element.getBBox();
            return {id, element, box, rank:meta.get(id)?.rank,
                horizontal_group:meta.get(id)?.horizontal_group, x:box.x+box.width/2,
                y:box.y+box.height/2, width:box.width};
        });
        const result = optimize(nodes,edges,clusters);
        svg.dataset.layoutMetrics=JSON.stringify(result);
        nodes.forEach(n=> {
            const dx=n.x-n.box.x-n.box.width/2;
            n.element.dataset.offsetX=dx;
            if(Number.isInteger(n.horizontal_group)) n.element.dataset.horizontalGroup=n.horizontal_group;
            n.element.setAttribute('transform', `translate(${dx},0)`);
        });
        // Refit the affine enclosures around the moved real nodes.
        clusters.forEach((cluster,i)=> {
            const group=[...svg.querySelectorAll('g.cluster')].find(g=>g.querySelector('title')?.textContent==='cluster_'+i);
            const members=nodes.filter(n=>cluster.nodes.includes(n.id));
            if (!group || !members.length) return;
            const x=Math.min(...members.map(n=>n.x-n.width/2))-12;
            const right=Math.max(...members.map(n=>n.x+n.width/2))+12;
            const top=Math.min(...members.map(n=>n.box.y))-30;
            const bottom=Math.max(...members.map(n=>n.box.y+n.box.height))+12;
            group.querySelectorAll('path,polygon').forEach(p=>p.remove());
            const rect=document.createElementNS('http://www.w3.org/2000/svg','rect');
            Object.entries({x,y:top,width:right-x,height:bottom-top,rx:10,fill:'none',stroke:cluster.color||'#6C5CE7','stroke-dasharray':'5,3'}).forEach(([k,v])=>rect.setAttribute(k,v));
            group.prepend(rect);
            group.querySelectorAll('text').forEach(t=> {t.setAttribute('x',(x+right)/2);t.setAttribute('y',top+15);});
        });
    }
    function reachable(seed, edges, reverse=false) {
        const found=new Set([seed]), queue=[seed], adjacency=new Map();
        edges.forEach(e=> {
            const backwards=reverse !== Boolean(e.reverseArrow);
            const a=backwards?e.target:e.source, b=backwards?e.source:e.target;
            if (!adjacency.has(a)) adjacency.set(a,[]);
            adjacency.get(a).push(b);
        });
        while(queue.length) for(const id of adjacency.get(queue.pop())||[]) {
            if (!found.has(id)) { found.add(id); queue.push(id); }
        }
        return found;
    }
    let focus = () => {};
    function fit(svg) {
        const graph=d3.select('#graph').graphviz(), bounds=[];
        svg.querySelectorAll('g.node, g.cluster, g.edge[id^="graph-edge-"]').forEach(g=> {
            const box=g.getBBox(), dx=Number(g.dataset.offsetX||0);
            bounds.push({x:box.x+dx,y:box.y,width:box.width,height:box.height});
        });
        if (!bounds.length) return;
        const left=Math.min(...bounds.map(b=>b.x))-25, right=Math.max(...bounds.map(b=>b.x+b.width))+25;
        const top=Math.min(...bounds.map(b=>b.y))-25, bottom=Math.max(...bounds.map(b=>b.y+b.height))+25;
        const width=svg.clientWidth, height=svg.clientHeight;
        const toolbar=document.querySelector('.graph-controls');
        const inset=(toolbar?.offsetHeight||40)+30;
        const scale=Math.min(width/(right-left),(height-inset)/(bottom-top));
        svg.setAttribute('viewBox',`0 0 ${width} ${height}`);
        graph.zoomSelection().call(graph.zoomBehavior().transform,
            d3.zoomIdentity.translate(width/2-scale*(left+right)/2,
                inset+(height-inset)/2-scale*(top+bottom)/2).scale(scale));
    }
    function interact(svg, metadata, edges) {
        const nodes=[...svg.querySelectorAll('g.node')];
        const mode=document.getElementById('graph-focus-mode'), status=document.getElementById('graph-focus-status');
        const witnessMode=document.getElementById('graph-witness-mode');
        focus = id => {
            const hierarchy=edges.filter(e=>e.hierarchy!==false);
            const up=reachable(id,hierarchy,true), down=reachable(id,hierarchy);
            const direct=new Set(edges.filter(e=>e.source===id||e.target===id).flatMap(e=>[e.source,e.target]));
            nodes.forEach(n=> {
                const key=n.querySelector('title').textContent;
                n.classList.toggle('graph-dimmed',!!id && !up.has(key) && !down.has(key) && !direct.has(key));
                n.classList.toggle('graph-selected',key===id);
            });
            edges.forEach((e,i)=>svg.querySelector('#graph-edge-'+i)?.classList.toggle('graph-dimmed',
                !!id && !(e.source===id || e.target===id || (up.has(e.source)&&up.has(e.target))||(down.has(e.source)&&down.has(e.target)))));
            status.textContent=id ? (metadata.find(n=>n.id===id)?.label||id).replaceAll('\\n',' ') : '';
        };
        document.getElementById('graph-clear-focus').onclick=()=>focus(null);
        document.getElementById('graph-fit').onclick=()=>fit(svg);
        mode.onchange=()=> {focus(null);status.textContent=mode.checked?'Click a node to focus; turn off to open cards.':'';};
        const graphGroup=svg.querySelector('g.graph');
        function updateLabels() {
            const matrix=graphGroup.getScreenCTM();
            const zoom=matrix ? Math.hypot(matrix.a,matrix.b) : 0;
            svg.classList.toggle('witness-labels-visible',witnessMode.value==='always'||(witnessMode.value==='auto'&&zoom>=0.85));
            svg.classList.toggle('witness-labels-hidden',witnessMode.value==='never');
        }
        witnessMode.onchange=updateLabels;
        new MutationObserver(updateLabels).observe(graphGroup,{attributes:true,attributeFilter:['transform']});
        new ResizeObserver(updateLabels).observe(svg);
        updateLabels();
        svg.querySelectorAll('g.edge').forEach(g=> {
            const path=g.querySelector('path');
            if (!path || !g.id.startsWith('graph-edge-')) return;
            const hit=path.cloneNode();
            hit.removeAttribute('marker-start');hit.removeAttribute('marker-end');
            hit.setAttribute('stroke','transparent');hit.setAttribute('stroke-width','12');
            hit.setAttribute('class','edge-hit-area');hit.setAttribute('pointer-events','stroke');
            g.insertBefore(hit,path);
        });
        // Curves and hit targets stay behind opaque parameter nodes.
        nodes.forEach(n=>graphGroup.appendChild(n));
        fit(svg);
        svg.dataset.readabilityReady='true';
    }
    return {metrics,optimize,place,interact,reachable,fit,focus:id=>focus(id)};
})();
