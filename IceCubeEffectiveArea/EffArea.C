void EffArea()
{
  TNtuple *n1 = new TNtuple("n1","n1","energy:exposure");
  n1->ReadFile("multi-year-diffuse-exposure-cosZenMin0.00-cosZenMax0.05.txt");

  TNtuple *n2 = new TNtuple("n2","n2","energy:exposure");
  n2->ReadFile("multi-year-diffuse-exposure-cosZenMin0.00-cosZenMax0.05-truncated-200TeV.txt");       

  n1->Draw("exposure:energy","","L");
  n2->Draw("exposure:energy","","L same");   

    
}
