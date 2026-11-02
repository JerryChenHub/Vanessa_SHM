using System;
using System.IO;
using System.Linq;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Runtime.InteropServices.ComTypes;
using SolidWorks.Interop.sldworks;

// Export workspace copies only. Do not save or rebuild source SLDPRT documents.
public static class ExportSkinGeometry {
 [DllImport("ole32.dll")] static extern int GetRunningObjectTable(int r,out IRunningObjectTable t);
 [DllImport("ole32.dll")] static extern int CreateBindCtx(int r,out IBindCtx c);
 static List<ISldWorks> Apps() {
  IRunningObjectTable t; IBindCtx c; GetRunningObjectTable(0,out t); CreateBindCtx(0,out c);
  IEnumMoniker e; t.EnumRunning(out e); var m=new IMoniker[1]; var apps=new List<ISldWorks>();
  while(e.Next(1,m,IntPtr.Zero)==0) { string n; m[0].GetDisplayName(c,null,out n);
   if(!n.StartsWith("SolidWorks_PID_"))continue;
   object v; t.GetObject(m[0],out v); apps.Add((ISldWorks)v);
  }
  return apps;
 }
 public static object Run(int pid,string folder) {
  var app=Apps().First(x=>x.GetProcessID()==pid);
  if(app.CommandInProgress)throw new Exception("SolidWorks is busy");
  var prior=(IModelDoc2)app.ActiveDoc; var results=new List<object>(); app.CommandInProgress=true;
  try {
   foreach(string tag in new[]{"1","2"}) {
    string path=Path.Combine(folder,tag+".SLDPRT"); IModelDoc2 doc=null; bool opened=false;
    try {
     if(app.GetOpenDocumentByName(path)!=null)throw new Exception("Export copy already open: "+path);
     int err=0,warn=0; doc=app.OpenDoc6(path,1,3,"",ref err,ref warn); opened=doc!=null;
     if(doc==null)throw new Exception("Open failed: "+err);
     int ae=0; app.ActivateDoc3(doc.GetTitle(),false,0,ref ae); doc.ClearSelection2(true);
     string output=Path.Combine(folder,tag+".x_t"); int se=0,sw=0;
     bool ok=doc.Extension.SaveAs(output,0,1,null,ref se,ref sw);
     if(!ok||!File.Exists(output))throw new Exception("Parasolid export failed: "+se);
     results.Add(new {part=tag,source=path,parasolid=output,open_error=err,open_warning=warn,save_error=se,save_warning=sw});
    } finally {if(opened&&doc!=null)app.CloseDoc(doc.GetTitle());}
   }
  } finally {if(prior!=null){int e=0;app.ActivateDoc3(prior.GetTitle(),false,0,ref e);} app.CommandInProgress=false;}
  return results;
 }
}
