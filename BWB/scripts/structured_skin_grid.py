"""CAD-projected, indexed span/chord grids for the supplied BWB geometry.

Topology is explicit: no paving, automatic node merge, or mesh smoothing.
Main-wing rows are constant CAD Z planes. Other patches use CAD-boundary
Coons interpolation followed by projection onto the actual trimmed surface.
The face/curve map is specific to the two checked-in Parasolid parts.
"""
from collections import defaultdict, deque
import numpy as np
from scipy.optimize import brentq
from femap_connection import check

# face, root chord, tip chord, leading edge, trailing-edge chain
# Signed curve IDs specify orientation from LE to TE or root to tip.
PATCHES = [
    ('main_lower', 2, [10], [4], [-7], [6,-26,-29], 'main'),
    ('main_upper', 3, [-8], [-3], [-7], [5,27,32], 'main'),
    ('transition_lower', 8, [4], [19], [15], [-31], 'transition'),
    ('transition_upper', 11, [-3], [-21], [15], [34], 'transition'),
    ('bend_lower', 7, [19], [18], [14], [-30,12], 'bend'),
    ('bend_upper', 10, [-21], [-20], [14], [33,28,17], 'bend'),
    ('winglet_lower', 6, [18], [23], [13], [11], 'winglet'),
    ('winglet_upper', 9, [-20], [-25], [13], [16], 'winglet'),
]

class Curves:
    def __init__(self, app):
        self.curve = app.feCurve
        self.cache = {}
        self.sampled = {}
        self.geometry_refs = {}
        self.endpoints = {}

    def point(self, cid, t):
        if cid < 0:
            cid, t = -cid, 1-t
        key = (cid, round(float(t),14))
        if key not in self.cache:
            check(self.curve.Get(cid), 'get CAD curve')
            rc, xyz = self.curve.ParamToXYZ(float(t))
            check(rc, 'evaluate CAD curve')
            self.cache[key] = np.asarray(xyz,dtype=float)
            position=tuple(np.round(self.cache[key]/1e-9).astype(np.int64))
            if t<1e-12 or t>1-1e-12:
                if cid not in self.endpoints:
                    rc,a,b=self.curve.EndPoints();check(rc,'curve endpoint IDs')
                    self.endpoints[cid]=(a,b)
                self.geometry_refs[position]=(3,self.endpoints[cid][0 if t<.5 else 1])
            elif self.geometry_refs.get(position,(5,0))[0]!=3:
                self.geometry_refs[position]=(4,cid)
        return self.cache[key].copy()

    def table(self, chain):
        key=tuple(chain)
        if key not in self.sampled:
            records=[]
            distance=0.
            for cid in chain:
                ts=np.linspace(0,1,129)
                points=np.array([self.point(cid,t) for t in ts])
                lengths=np.r_[0,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))]
                records.append((cid,ts,lengths,points,distance))
                distance+=lengths[-1]
            self.sampled[key]=(records,distance)
        return self.sampled[key]

    def at(self, chain, f, coordinate=None):
        if f <= 1e-14:return self.point(chain[0],0)
        if f >= 1-1e-14:return self.point(chain[-1],1)
        records,length=self.table(chain)
        if coordinate is None:
            wanted=f*length
            for cid,ts,ls,points,start in records:
                if wanted<=start+ls[-1]+1e-13:
                    return self.point(cid,np.interp(wanted-start,ls,ts))
        else:
            start=self.point(chain[0],0)[coordinate]
            end=self.point(chain[-1],1)[coordinate]
            target=(1-f)*start+f*end
            for cid,ts,ls,points,d in records:
                vals=points[:,coordinate]-target
                crossings=np.flatnonzero(vals[:-1]*vals[1:]<=0)
                if len(crossings):
                    k=int(crossings[0])
                    t=brentq(lambda v:self.point(cid,v)[coordinate]-target,ts[k],ts[k+1],xtol=1e-13)
                    return self.point(cid,t)
        raise ValueError(f'Cannot sample boundary {chain}, fraction={f}')

def _orient_consistently(nodes, elements):
    """Orient connected shells by shared edges, then outward using signed volume."""
    uses=defaultdict(list)
    for i,e in enumerate(elements):
        n=e['nodes']
        for a,b in zip(n,n[1:]+n[:1]):
            uses[tuple(sorted((a,b)))].append((i,1 if a<b else -1))
    if any(len(u)>2 for u in uses.values()):
        raise RuntimeError('Structured grid has a nonmanifold edge')
    neighbours=defaultdict(list)
    for u in uses.values():
        if len(u)==2:
            (a,da),(b,db)=u
            neighbours[a].append((b,-da*db));neighbours[b].append((a,-da*db))
    sign={0:1};queue=deque([0])
    while queue:
        a=queue.popleft()
        for b,relation in neighbours[a]:
            wanted=sign[a]*relation
            if b in sign:
                if sign[b]!=wanted:raise RuntimeError('Inconsistent surface orientation')
            else:sign[b]=wanted;queue.append(b)
    if len(sign)!=len(elements):raise RuntimeError('Disconnected structured patches')
    for i,e in enumerate(elements):
        if sign[i]<0:e['nodes'].reverse()
    # Root lies in a plane through the first grid node. A virtual planar root cap
    # therefore contributes zero to the signed-volume sum about this origin.
    # It is used ONLY to orient normals; no root cap elements are generated.
    centre=nodes[min(nodes)]
    volume=0.
    for e in elements:
        p=np.array([nodes[n] for n in e['nodes']])-centre
        for j in range(1,len(p)-1):volume+=np.dot(p[0],np.cross(p[j],p[j+1]))/6
    if volume<0:
        for e in elements:e['nodes'].reverse()
    return abs(float(volume))

def generate(app, size=.020, chord_divisions=24):
    curves=Curves(app);surface=app.feSurface
    nodes={};node_lookup={};elements=[];grids={};stats={};node_faces=defaultdict(set)
    def node(xyz,sid):
        # Only identical CAD boundary samples merge. Resolution is 1 nanometre.
        key=tuple(np.round(np.asarray(xyz)/1e-9).astype(np.int64))
        nid=node_lookup.get(key)
        if nid is None:
            nid=len(nodes)+1;nodes[nid]=np.asarray(xyz,dtype=float);node_lookup[key]=nid
        elif np.linalg.norm(nodes[nid]-xyz)>2e-9:
            raise RuntimeError('Boundary node identity tolerance exceeded')
        node_faces[nid].add(sid)
        return nid
    def element(conn,sid,patch):
        clean=[]
        for n in conn:
            if n not in clean:clean.append(int(n))
        if len(clean)<3:raise RuntimeError('Collapsed grid element')
        elements.append(dict(nodes=clean,surface_id=sid,patch=patch))

    sections={'main':max(2,round(.78/size)), 'transition':max(3,round(.035/size)),
              'bend':max(8,round(.065/size)), 'winglet':max(2,round(.185/size))}
    for name,sid,root,tip,le,te,zone in PATCHES:
        check(surface.Get(sid),f'load {name}')
        ns=sections[zone];nc=chord_divisions
        R=np.array([curves.at(root,f,0) for f in np.linspace(0,1,nc+1)])
        T=np.array([curves.at(tip,f,0) for f in np.linspace(0,1,nc+1)])
        L=np.array([curves.at(le,f,2 if zone=='main' else None) for f in np.linspace(0,1,ns+1)])
        E=np.array([curves.at(te,f,2 if zone=='main' else None) for f in np.linspace(0,1,ns+1)])
        grid=np.zeros((ns+1,nc+1),dtype=int)
        projected_count=0;max_projection=0.;inside_failures=[]
        for j in range(ns+1):
            v=j/ns
            span=(L[min(j+1,ns)]+E[min(j+1,ns)]-L[max(j-1,0)]-E[max(j-1,0)])/2
            normal=np.cross(span,E[j]-L[j]);normal/=np.linalg.norm(normal)
            if zone=='main':normal=np.array([0.,1.,0.])
            for i in range(nc+1):
                u=i/nc
                if j==0:p=R[i]
                elif j==ns:p=T[i]
                elif i==0:p=L[j]
                elif i==nc:p=E[j]
                else:
                    corners=(1-u)*(1-v)*R[0]+u*(1-v)*R[-1]+(1-u)*v*T[0]+u*v*T[-1]
                    seed=(1-v)*R[i]+v*T[i]+(1-u)*L[j]+u*E[j]-corners
                    if zone=='main':
                        seed[0]=(1-u)*L[j,0]+u*E[j,0]
                        seed[2]=L[j,2]
                    rc,xyz=surface.ProjectAlongVector(True,tuple(normal),tuple(seed))
                    check(rc,f'project {name} row {j} col {i}')
                    p=np.asarray(xyz,dtype=float)
                    if surface.InsideXYZ(tuple(p))!=-1:
                        inside_failures.append([j,i])
                    if zone=='main' and max(abs(p[[0,2]]-seed[[0,2]]))>1e-8:
                        raise RuntimeError('Main-wing projection moved the parallel station')
                    max_projection=max(max_projection,float(np.linalg.norm(p-seed)))
                    projected_count+=1
                grid[j,i]=node(p,sid)
        if inside_failures:
            raise RuntimeError(f'{name}: projected nodes outside trimmed face: {inside_failures[:10]}')
        for j in range(ns):
            for i in range(nc):element([grid[j,i],grid[j,i+1],grid[j+1,i+1],grid[j+1,i]],sid,name)
        grids[name]=grid.tolist()
        stats[name]=dict(surface_id=sid,span_divisions=ns,chord_divisions=nc,
                         projected_nodes=projected_count,max_projection_m=max_projection,
                         main_row_z_range_max_m=max(float(np.ptp([nodes[n][2] for n in row])) for row in grid) if zone=='main' else None)
        print(f'  {name}: {ns} x {nc} quads; projection max {max_projection:.3g} m',flush=True)

    # Close the trailing edge with matching upper/lower row nodes. The tiny CAD
    # face boundaries are not treated as mesh boundaries; no new internal ribs.
    te_faces={'main':[4,16,15], 'transition':[15], 'bend':[15,13,14], 'winglet':[12]}
    for zone in sections:
        low=np.array(grids[f'{zone}_lower']);high=np.array(grids[f'{zone}_upper'])
        for j in range(len(low)-1):
            conn=[low[j,-1],high[j,-1],high[j+1,-1],low[j+1,-1]]
            centre=np.mean([nodes[int(n)] for n in conn],axis=0)
            closest=[]
            for sid in te_faces[zone]:
                surface.Get(sid);rc,p=surface.ProjectToClosest(tuple(centre))
                if rc==-1:closest.append((np.linalg.norm(np.asarray(p)-centre),sid))
            sid=min(closest)[1]
            element(conn,sid,f'{zone}_trailing_edge')
    # External winglet tip closure: one unavoidable triangle where LE thickness
    # tends to zero, then an ordered strip of quads across the chord.
    low=np.array(grids['winglet_lower'])[-1];high=np.array(grids['winglet_upper'])[-1]
    for i in range(chord_divisions):element([low[i],low[i+1],high[i+1],high[i]],5,'tip_cap')
    volume=_orient_consistently(nodes,elements)

    # Explicit row/column topology is retained for QA and subsequent rib insertion.
    node_geometry={str(nid):curves.geometry_refs.get(tuple(np.round(p/1e-9).astype(np.int64)),
                                                    (5,min(node_faces[nid]))) for nid,p in nodes.items()}
    report=dict(method='CAD-projected structured span/chord grids',chord_divisions=chord_divisions,
                span_divisions=sections,patches=stats,patch_node_ids=grids,
                virtual_root_cap_volume_m3=volume,root_cap_meshed=False,
                node_geometry=node_geometry,
                trailing_edge_representation='ruled upper/lower edge strips; tiny CAD seams simplified')
    return nodes,elements,report

def put_in_femap(app,nodes,elements,node_geometry):
    ordered=sorted(nodes)
    check(app.feNode.PutCoordArray(len(ordered),tuple(ordered),tuple(float(v) for n in ordered for v in nodes[n])), 'create structured nodes')
    check(app.feMeshAttachNodes(len(ordered),tuple(ordered),
                               tuple(node_geometry[str(n)][0] for n in ordered),
                               tuple(node_geometry[str(n)][1] for n in ordered)), 'associate structured nodes to CAD')
    e=app.feElem
    grouped=defaultdict(list)
    for eid,entry in enumerate(elements,1):
        e.type=17;e.topology=4 if len(entry['nodes'])==4 else 2;e.propID=1
        e.vnode=tuple(entry['nodes'])+(0,)*(20-len(entry['nodes']))
        e.geomID=entry['surface_id']
        check(e.Put(eid),'create structured shell')
        entry['id']=eid
        grouped[entry['surface_id']].append(eid)
    for sid,eids in grouped.items():
        sel=app.feSet;sel.AddArray(len(eids),tuple(eids))
        check(app.feMeshAttach(8,sel.ID,5,sid,False),'associate shells to CAD')
    return dict(grouped)
