"use strict";
const fs = require("node:fs");
const path = require("node:path");
const crypto = require("node:crypto");
const root = process.cwd();
const fixtureDir = path.join(root, "pgfrag/experiments/explicit_fragment_fixture_20261006");
const outDir = path.join(root, "pgfrag/tmp/review_20261006/explicit_fragment_fixture_20261006");
const auditDir = path.join(root, "pgfrag/tmp/review_20261006/fragment_geometry_export_independent");
const input = JSON.parse(fs.readFileSync(path.join(fixtureDir, "input_fixture.json"), "utf8"));
const data = JSON.parse(fs.readFileSync(path.join(outDir, "explicit_fragment_fixture_output.json"), "utf8"));
const ledger = JSON.parse(fs.readFileSync(path.join(outDir, "interface_pair_ledger.json"), "utf8"));
const savedChecks = JSON.parse(fs.readFileSync(path.join(outDir, "checks.json"), "utf8"));
const add=(a,b)=>a.map((v,i)=>v+b[i]), sub=(a,b)=>a.map((v,i)=>v-b[i]);
const mul=(a,s)=>a.map(v=>v*s), dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const norm=a=>Math.sqrt(dot(a,a)), area=(a,b,c)=>norm(cross(sub(b,a),sub(c,a)))/2;
const unit=a=>mul(a,1/norm(a));
const coordKey=p=>p.map(v=>Math.round(v*1e9)).join(",");
const nearly=(a,b,t=1e-9)=>Math.abs(a-b)<=t*Math.max(1,Math.abs(a),Math.abs(b));
const assert=(ok,msg)=>{if(!ok) throw Error(msg);};
const hashFile=p=>crypto.createHash("sha256").update(fs.readFileSync(p)).digest("hex");
const sourceFiles=[
 "pgfrag/experiments/explicit_fragment_fixture_20261006/input_fixture.json",
 "pgfrag/experiments/explicit_fragment_fixture_20261006/explicit_fragment_fixture.js",
 "pgfrag/tmp/review_20261006/explicit_fragment_fixture_20261006/explicit_fragment_fixture_output.json",
 "pgfrag/tmp/review_20261006/explicit_fragment_fixture_20261006/interface_pair_ledger.json",
 "pgfrag/tmp/review_20261006/explicit_fragment_fixture_20261006/checks.json",
 "pgfrag/tmp/review_20261006/explicit_fragment_fixture_20261006/piece_negative.obj",
 "pgfrag/tmp/review_20261006/explicit_fragment_fixture_20261006/piece_positive.obj"
];
function objRead(side) {
 const file=path.join(outDir,"piece_"+side+".obj"), vertices=[], faces=[], groups=[];
 let active="";
 for(const line of fs.readFileSync(file,"utf8").split(/\r?\n/)) {
  const t=line.trim().split(/\s+/);
  if(t[0]==="v") vertices.push(t.slice(1,4).map(Number));
  else if(t[0]==="g") active=t[1]||"";
  else if(t[0]==="f") { faces.push(t.slice(1).map(x=>Number(x.split("/")[0])-1)); groups.push(active); }
 }
 return {vertices,faces,groups};
}
function inspectMesh(mesh) {
 const edgeMap=new Map(), adj=Array.from({length:mesh.vertices_m.length},()=>new Set());
 let signedVolume=0, moment=[0,0,0], surfaceArea=0;
 const lineageAreas={};
 for(let fi=0;fi<mesh.triangles.length;fi++) {
  const tri=mesh.triangles[fi], pts=tri.map(i=>mesh.vertices_m[i]);
  assert(tri.length===3&&new Set(tri).size===3,"invalid triangle indices");
  const a=pts[0],b=pts[1],c=pts[2], ar=area(a,b,c);
  assert(ar>0,"degenerate output triangle");
  surfaceArea+=ar;
  const tv=dot(a,cross(b,c))/6;
  signedVolume+=tv; moment=add(moment,mul(add(add(a,b),c),tv/4));
  const lin=mesh.face_lineage[fi];
  const tag=lin.kind==="original_surface_face"?lin.original_surface_face_id:lin.surface_id;
  lineageAreas[tag]=(lineageAreas[tag]||0)+ar;
  for(let k=0;k<3;k++) {
   const x=tri[k],y=tri[(k+1)%3],key=x<y?x+","+y:y+","+x;
   if(!edgeMap.has(key))edgeMap.set(key,[]);
   edgeMap.get(key).push([x,y]); adj[x].add(y); adj[y].add(x);
  }
 }
 const badEdges=[...edgeMap].filter(([,uses])=>uses.length!==2||uses[0][0]!==uses[1][1]||uses[0][1]!==uses[1][0]);
 const seen=new Set(mesh.vertices_m.length?[0]:[]), stack=[...seen];
 while(stack.length){for(const j of adj[stack.pop()])if(!seen.has(j)){seen.add(j);stack.push(j)}}
 const euler=mesh.vertices_m.length-edgeMap.size+mesh.triangles.length;
 assert(!badEdges.length,mesh.side+": bad directed edge incidence ("+badEdges.length+")");
 assert(seen.size===mesh.vertices_m.length,mesh.side+": disconnected or isolated vertices");
 assert(euler===2,mesh.side+": Euler "+euler);
 assert(signedVolume>0,mesh.side+": signed volume is not positive");
 const centroid=mul(moment,1/signedVolume);
 assert(nearly(signedVolume,mesh.metrics.volume_m3,1e-12),mesh.side+": saved volume mismatch");
 for(let k=0;k<3;k++)assert(nearly(centroid[k],mesh.metrics.centroid_m[k],1e-11),mesh.side+": saved centroid mismatch");
 const obj=objRead(mesh.side);
 assert(obj.vertices.length===mesh.vertices_m.length&&obj.faces.length===mesh.triangles.length,mesh.side+": OBJ count mismatch");
 for(let i=0;i<obj.vertices.length;i++)for(let k=0;k<3;k++)assert(nearly(obj.vertices[i][k],mesh.vertices_m[i][k],1e-15),mesh.side+": OBJ vertex mismatch");
 const normGroup=s=>s.replace(/[^A-Za-z0-9_.-]/g,"_");
 for(let i=0;i<obj.faces.length;i++) {
  assert(obj.faces[i].join(",")===mesh.triangles[i].join(","),mesh.side+": OBJ face "+i+" differs from JSON");
  const lin=mesh.face_lineage[i];
  const expected=normGroup(lin.kind==="original_surface_face"?lin.original_surface_face_id:lin.surface_id+"_fixture_cut");
  assert(obj.groups[i]===expected,mesh.side+": OBJ group/lineage mismatch at face "+i);
 }
 return {side:mesh.side,vertex_count:mesh.vertices_m.length,face_count:mesh.triangles.length,edge_count:edgeMap.size,euler,single_connected_component:true,closed_consistently_oriented:true,signed_volume_m3:signedVolume,saved_volume_m3:mesh.metrics.volume_m3,centroid_m:centroid,saved_centroid_m:mesh.metrics.centroid_m,surface_area_m2:surfaceArea,lineage_area_m2:lineageAreas,obj_matches_json_geometry_and_face_groups:true,obj_sha256:hashFile(path.join(outDir,"piece_"+mesh.side+".obj"))};
}
function polygonArea(poly) {
 if(poly.length<3)return 0;
 let v=[0,0,0]; for(let i=0;i<poly.length;i++)v=add(v,cross(poly[i],poly[(i+1)%poly.length]));
 return norm(v)/2;
}
function clip(poly,n,d,keep) {
 let out=[];
 const s=p=>dot(n,p)-d, inside=p=>keep==="negative"?s(p)<=1e-9:s(p)>=-1e-9;
 let prev=poly[poly.length-1], pin=inside(prev);
 for(const cur of poly) {
  const cin=inside(cur);
  if(cin!==pin){const dp=s(prev),dc=s(cur),t=dp/(dp-dc);out.push(add(prev,mul(sub(cur,prev),t)))}
  if(cin)out.push(cur);
  prev=cur;pin=cin;
 }
 const clean=[];
 for(const p of out)if(!clean.length||norm(sub(clean[clean.length-1],p))>1e-9)clean.push(p);
 if(clean.length>1&&norm(sub(clean[0],clean[clean.length-1]))<=1e-9)clean.pop();
 return clean;
}
const meshes=Object.fromEntries(data.pieces.map(m=>[m.side,inspectMesh(m)]));
const n=data.cut_surface.plane.normal,d=data.cut_surface.plane.offset_m;
const body=input.original_body, surf=input.cut_surface;
const sign=body.vertices_m.map(p=>dot(n,p)-d);
assert(Math.min(...sign)<-1e-6&&Math.max(...sign)>1e-6,"input plane does not strictly split original vertices");
const edgeSet=new Map();
for(const f of body.faces)for(let i=0;i<f.vertices.length;i++){
 const a=f.vertices[i],b=f.vertices[(i+1)%f.vertices.length],k=a<b?a+","+b:b+","+a;
 if(!edgeSet.has(k))edgeSet.set(k,[a,b]);
}
const hits=[];
for(const [a,b] of edgeSet.values())if((sign[a]<0&&sign[b]>0)||(sign[a]>0&&sign[b]<0)){
 const t=sign[a]/(sign[a]-sign[b]);hits.push(add(body.vertices_m[a],mul(sub(body.vertices_m[b],body.vertices_m[a]),t)));
}
const hitByKey=new Map(hits.map(p=>[coordKey(p),p]));
const declaredByKey=new Map(surf.vertices_m.map(p=>[coordKey(p),p]));
assert(hitByKey.size===declaredByKey.size&&[...hitByKey.keys()].every(k=>declaredByKey.has(k)),"cut patch vertices do not equal independently reconstructed body-edge intersections");
assert(surf.vertices_m.every(p=>Math.abs(dot(n,p)-d)<1e-9),"cut vertices not on declared plane");
const inputTriArea=surf.triangles.reduce((sum,t)=>sum+area(...t.map(i=>surf.vertices_m[i])),0);
assert(nearly(inputTriArea,data.cut_surface.area_m2,1e-12),"declared cut-surface area mismatch");
const expectedCutBoundary=[];
const triEdges=new Map();
for(const t of surf.triangles)for(let i=0;i<3;i++){
 const a=t[i],b=t[(i+1)%3],k=a<b?a+","+b:b+","+a;
 if(!triEdges.has(k))triEdges.set(k,[]);
 triEdges.get(k).push([a,b]);
}
for(const uses of triEdges.values())if(uses.length===1)expectedCutBoundary.push(uses[0]);
assert(expectedCutBoundary.length===surf.vertices_m.length,"cut patch boundary is not one simple cycle");
for(const [a,b] of expectedCutBoundary) {
 const pa=surf.vertices_m[a],pb=surf.vertices_m[b];
 assert([0,1,2].some(k=>Math.abs(pa[k]-pb[k])<1e-9&&Math.abs(Math.abs(pa[k])-1)<1e-9),"cut boundary edge is not on the input cube exterior");
}
const sideLineage={negative:{},positive:{}};
for(const mesh of data.pieces) {
 const mside=mesh.side;
 for(let i=0;i<mesh.face_lineage.length;i++){
  const lin=mesh.face_lineage[i], tri=mesh.triangles[i].map(v=>mesh.vertices_m[v]);
  if(lin.kind==="original_surface_face"){
   const f=body.faces.find(x=>x.id===lin.original_surface_face_id);
   assert(!!f,"unknown parent face lineage "+lin.original_surface_face_id);
   const src=f.vertices.map(v=>body.vertices_m[v]);
   let axis=-1; for(let k=0;k<3;k++)if(src.every(p=>Math.abs(p[k]-src[0][k])<1e-10))axis=k;
   assert(axis>=0,"fixture parent face "+f.id+" not axis planar");
   for(const p of tri) {
    assert(Math.abs(p[axis]-src[0][axis])<1e-9,"face "+i+" leaves claimed parent face plane");
    for(let k=0;k<3;k++)assert(p[k]>=Math.min(...src.map(q=>q[k]))-1e-9&&p[k]<=Math.max(...src.map(q=>q[k]))+1e-9,"face "+i+" leaves claimed parent face bounds");
   }
   sideLineage[mside][lin.original_surface_face_id]=(sideLineage[mside][lin.original_surface_face_id]||0)+area(...tri);
  } else {
   assert(lin.kind==="explicit_fracture_face_fixture"&&lin.surface_id===surf.surface_id&&lin.fixture_geometry_only===true,"invalid cut face lineage at "+mside+":"+i);
   assert(Number.isInteger(lin.input_surface_triangle_index)&&surf.triangles[lin.input_surface_triangle_index],"invalid input cut triangle lineage");
  }
 }
 const expectedOriginalFaces=Object.fromEntries(body.faces.map(f=>[f.id,polygonArea(clip(f.vertices.map(v=>body.vertices_m[v]),n,d,mside))]));
 for(const f of body.faces)assert(nearly(sideLineage[mside][f.id]||0,expectedOriginalFaces[f.id],1e-10),mside+": clipped source-face area mismatch for "+f.id);
}
const pairChecks=[];
for(const pair of ledger.triangle_pairs) {
 const neg=data.pieces.find(m=>m.side==="negative"),pos=data.pieces.find(m=>m.side==="positive");
 const nt=neg.triangles[pair.negative_face_index],pt=pos.triangles[pair.positive_face_index];
 assert(pair.negative_face_id===neg.piece_id+":f"+pair.negative_face_index&&pair.positive_face_id===pos.piece_id+":f"+pair.positive_face_index,"ledger face IDs/index disagree");
 const npts=nt.map(i=>neg.vertices_m[i]),ppts=pt.map(i=>pos.vertices_m[i]);
 const keySet=ps=>ps.map(coordKey).sort().join("|");
 assert(keySet(npts)===keySet(ppts)&&keySet(npts)===keySet(pair.mother_coordinates_m),"paired cut faces do not share actual exported mother coordinates");
 const nn=unit(cross(sub(npts[1],npts[0]),sub(npts[2],npts[0])));
 const pn=unit(cross(sub(ppts[1],ppts[0]),sub(ppts[2],ppts[0])));
 const na=area(...npts),pa=area(...ppts),la=area(...pair.mother_coordinates_m);
 assert(dot(nn,pair.negative_outward_normal)>1-1e-12&&dot(pn,pair.positive_outward_normal)>1-1e-12,"ledger normal does not match exported triangle");
 assert(dot(nn,pn)<-1+1e-12&&dot(nn,n)>1-1e-12&&dot(pn,n)<-1+1e-12,"actual cap normals not oppositely outward");
 assert(nearly(na,pa,1e-12)&&nearly(na,la,1e-12),"triangle pair area differs");
 pairChecks.push({pair_id:pair.pair_id,negative_area_m2:na,positive_area_m2:pa,normal_dot:dot(nn,pn),coords_match:true});
}
for(const side of ["negative","positive"]){
 const mesh=data.pieces.find(m=>m.side===side), cutLin=mesh.face_lineage.filter(l=>l.kind==="explicit_fracture_face_fixture");
 const inds=cutLin.map(l=>l.input_surface_triangle_index).sort((a,b)=>a-b);
 assert(inds.join(",")===surf.triangles.map((_,i)=>i).join(","),side+": input cut triangle lineage incomplete/duplicated");
 const cutArea=meshes[side].lineage_area_m2[surf.surface_id]||0;
 assert(nearly(cutArea,inputTriArea,1e-12),side+": cut face area not equal input surface");
 const maxSigned=Math.max(...mesh.vertices_m.map(p=>dot(n,p)-d));
 const minSigned=Math.min(...mesh.vertices_m.map(p=>dot(n,p)-d));
 assert(side==="negative"?maxSigned<=1e-8:minSigned>=-1e-8,side+": vertices cross cut halfspace");
}
const outputVolumes=meshes.negative.signed_volume_m3+meshes.positive.signed_volume_m3;
const originalVolume=8;
assert(nearly(outputVolumes,originalVolume,1e-12),"output pieces do not conserve original cube volume");
assert(nearly(data.original_body.volume_m3,originalVolume,1e-12),"saved original volume differs from cube volume");
assert(nearly(data.pieces.reduce((s,m)=>s+m.metrics.mass_kg,0),originalVolume*body.density_kg_m3,1e-12),"piece mass does not conserve parent mass");
const expectedExternalArea=body.faces.reduce((sum,f)=>sum+polygonArea(f.vertices.map(v=>body.vertices_m[v])),0);
const actualExternalArea=Object.values(meshes).reduce((s,m)=>s+Object.entries(m.lineage_area_m2).filter(([k])=>k.startsWith("original:")).reduce((z,[,v])=>z+v,0),0);
assert(nearly(actualExternalArea,expectedExternalArea,1e-10),"original-face lineage does not conserve parent exterior area");
const receipt={
 schema:"pgfrag.explicit_fragment_fixture_independent_readback.v1",
 generated_utc:new Date().toISOString(),
 method:"One bounded independent Node.js readback; no solver or experiment test executed.",
 input_qualification:data.qualification,mechanics_qualification:data.mechanics_qualification,physical_fracture_claim:data.physical_fracture_claim,
 source_sha256:Object.fromEntries(sourceFiles.map(rel=>[rel,hashFile(path.join(root,rel))])),
 geometry:{original_body_id:body.body_id,original_vertex_count:body.vertices_m.length,original_face_count:body.faces.length,cut_surface_id:surf.surface_id,cut_surface_kind:surf.geometry_kind,cut_surface_vertex_count:surf.vertices_m.length,cut_surface_triangle_count:surf.triangles.length,reconstructed_original_edge_intersections_m:[...hitByKey.values()],input_cut_area_m2:inputTriArea,independently_reconstructed_full_section:true,cut_patch_boundary_on_cube_exterior:true,piece_meshes:meshes,actual_interface_triangle_pairs:pairChecks,paired_interface_area_negative_m2:meshes.negative.lineage_area_m2[surf.surface_id],paired_interface_area_positive_m2:meshes.positive.lineage_area_m2[surf.surface_id],original_exterior_area_m2:expectedExternalArea,lineaged_original_exterior_area_m2:actualExternalArea,independent_piece_volume_sum_m3:outputVolumes,original_box_volume_m3:originalVolume,mass_sum_kg:data.pieces.reduce((s,m)=>s+m.metrics.mass_kg,0),area_lineage_matches_clipped_source_faces:true,conservation_and_halfspace_checks:true},
 saved_producer_checks:savedChecks,
 findings:[
  "Geometry, orientation, closure, OBJ/JSON agreement, paired coordinates/normals/areas, and source-face/cut-face labels independently agree for this supplied planar cube fixture.",
  "These checks validate a deterministic explicit planar split; they do not identify a crack from broken bonds or establish material-realistic rough fracture geometry."
 ]
};
fs.writeFileSync(path.join(auditDir,"independent_readback.json"),JSON.stringify(receipt,null,2)+"\n","utf8");
process.stdout.write(JSON.stringify({path:path.join(auditDir,"independent_readback.json"),output_volumes_m3:[meshes.negative.signed_volume_m3,meshes.positive.signed_volume_m3],cut_area_m2:inputTriArea,external_area_m2:actualExternalArea,readback:"passed"},null,2)+"\n");

