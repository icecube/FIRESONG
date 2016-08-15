/* The following flux corresponds to
E^2dN/dE = 10^-8 (E/100 TeV)^-0.1 [GeV/cm^2.s.sr]
when assigned to a *single* source in the sky

 */

Double_t flux(Double_t energy)
{
  Double_t A = 3.9738353e-3; // Units 1/(GeV.m^2.s)
  Double_t alpha = -2.1;
  return A*TMath::Power(energy,alpha);
}

/* 
multi-year-diffuse-exposure-cosZenMin0.05-cosZenMax0.10-truncated-200TeV.txt 84.261 87.134 105.896
multi-year-diffuse-exposure-cosZenMin0.00-cosZenMax0.05-truncated-200TeV.txt 87.134 90.00  185.11
multi-year-diffuse-exposure-cosZenMin-0.05-cosZenMax0.00-truncated-200TeV.txt 90.00 92.866 177.775
multi-year-diffuse-exposure-cosZenMin-0.10-cosZenMax-0.05-truncated-200TeV.txt 92.866 95.739 155.445
multi-year-diffuse-exposure-cosZenMin-0.15-cosZenMax-0.10-truncated-200TeV.txt 95.739 98.627 133.492
multi-year-diffuse-exposure-cosZenMin-0.20-cosZenMax-0.15-truncated-200TeV.txt 98.627 101.537 113.77
multi-year-diffuse-exposure-cosZenMin-0.25-cosZenMax-0.20-truncated-200TeV.txt 101.537 104.478 99.8513
multi-year-diffuse-exposure-cosZenMin-0.30-cosZenMax-0.25-truncated-200TeV.txt 104.478 107.458 87.4213
multi-year-diffuse-exposure-cosZenMin-0.35-cosZenMax-0.30-truncated-200TeV.txt 107.458 110.487 76.7264
multi-year-diffuse-exposure-cosZenMin-0.40-cosZenMax-0.35-truncated-200TeV.txt 110.487 113.578 63.8566
*/

void SignalExpectation()
{
  TNtuple *n = new TNtuple("n","n","e:a");
  n->ReadFile("multi-year-diffuse-exposure-cosZenMin0.00-cosZenMax0.05.txt");

  Float_t energy;
  Float_t exposure;
  Float_t nevents=0;
  Float_t nevents_total = 0;
  n->SetBranchAddress("e",&energy);
  n->SetBranchAddress("a",&exposure);

  for (int i=0;i<n->GetEntries();i++)
    {
      n->GetEntry(i);      
      Double_t L10E_low = log10(energy);
      Double_t L10E_hi = L10E_low + 0.2;
      Double_t L10E_center = L10E_low + 0.1;
      Double_t energy_center = TMath::Power(10,L10E_center);
      Double_t energy_hi = TMath::Power(10,L10E_hi);
      nevents = flux(energy_center) * exposure * (energy_hi-energy);
      nevents_total += nevents;
      cout << i << " " 
	   << energy << " "
	   << energy_center <<  " "
	   << energy_hi << " "
	   << exposure << " "
	   << flux(energy_center) << " "
	   << nevents << endl;
    }
  cout << nevents_total << endl;
}
