{
  // D0 toy mc
  //double x[10]  = {1.05   , 1.15   , 1.25   , 1.35   , 1.45   , 1.55   , 1.65   , 1.75   , 1.85,  , 1.95};
  //double y[10]  = {1.86356, 1.86352, 1.86382, 1.86416, 1.86461, 1.86498, 1.86547, 1.86586, 1.86644, 1.86720};
  //double y2[10] = {1.86344, 1.86370, 1.86396, 1.86421, 1.86456, 1.86500, 1.86546, 1.86599, 1.86677, 1.86724};
  //double xe[10] = {0};
  //double ye[10] = {6.30e-4, 2.13e-4, 1.15e-4, 7.4e-5 , 6.0e-5 , 6.1e-5 , 7.5e-5 , 1.14e-4, 2.19e-4, 5.84e-4};
  

  // Ks toy mc
  double x[9]  = {    0.55,     0.65,     0.75,     0.85,     0.95,      1.1,      1.3,      1.5,      1.7};
  double y2[9] = {0.497483, 0.497440, 0.497497, 0.497506, 0.497532, 0.497596, 0.497701, 0.497829, 0.498041};
  double y[9]  = {0.497446, 0.497591, 0.497479, 0.497478, 0.497527, 0.497594, 0.497708, 0.497799, 0.498109};
  double ye[9] = {1.08e-4,  5.9e-5,   3.58e-5,  2.55e-5,  1.99e-5,  1.26e-5,  1.56e-5,  3.16e-5,  9.53e-5};
  double xe[9] = {0};

  TGraphErrors *graph = new TGraphErrors(9,x,y,xe,ye);
  graph->SetLineColor(2);
  graph->SetMarkerStyle(5);
  graph->SetMarkerColor(2);
  graph->SetFillColor(0);
  graph->SetTitle("Compare mass spectrum fit peak with estimate from formula");
  graph->GetXaxis()->SetTitle("p1(GeV/c)");
  graph->GetYaxis()->SetTitle("mass(GeV/c^{2})");

  TGraph *graph2= new TGraph(9,x,y2);
  graph2->SetLineColor(3);
  graph2->SetMarkerStyle(5);
  graph2->SetMarkerColor(3);
  graph2->SetFillColor(0);
  
  graph->Draw("AP");
  graph2->Draw("P");

  TLegend *legend = new TLegend(0.15,0.75,0.45,0.85);
  legend->AddEntry(graph,"mass spectrum fit result");
  legend->AddEntry(graph2,"mass peak estimate");
  legend->Draw();
}
