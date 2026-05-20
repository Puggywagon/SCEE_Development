Notes:
There are big deviations between the results I have obtained for Ethanol using OPLS-AA, TraPPE-UA and Cecilias PolCA. Now I dunno if this is specifically a AA thing, an opls thing or what. What I am going to do is run test simulations comparing different AA and UA models.



**Preparation**
- [ ] Look through Leo's Script folder
	- [x] Refresh on how the new scripts work compared to the old scripts ✅ 2026-05-20
	- [x] Make changes to new scripts based on changes made to the old scripts ✅ 2026-05-20
	- [x] Make Changes to Scripts ✅ 2026-05-20
		- [x] To allow to work with AA and UA models (i.e. the central molecule is dif from the solvent) ✅ 2026-05-20
		- [x] Make changes to where the script runs ✅ 2026-05-20
			- [x] All scripts stored in a home base ✅ 2026-05-20
			- [x] Molecule of interest in a specified directory where all outputs for that molecule will be saved ✅ 2026-05-20
			- [x] Save figures, csvs and output data into a data folder. ✅ 2026-05-20
		- [x] Changes based on the previous scripts changes ✅ 2026-05-20
			- [x] Basis Sets for all gaussian calcs ✅ 2026-05-20
			- [x] Adjust Scale ratio to have a low dip, dip and high dip (0.8,1,1.5 rather than 1,1.35,1.7) ✅ 2026-05-20
			- [x] Allow for more replicas ✅ 2026-05-20
		- [ ] Analysis additions
			- [x] Make a step for generating the ndxs of the molecule (going to not do this) ✅ 2026-05-20
			- [x] Perform the RDFs for the molecule (going to not do this) ✅ 2026-05-20
			- [x] Calculate the coordination number (going to not do this)  ✅ 2026-05-20
			- [x] calculate the enthalpy of vapourisation ✅ 2026-05-20
			- [x] Calculate the heat capacity (going to not do this)  ✅ 2026-05-20
			- [x] Generate associated figures ✅ 2026-05-20
- [ ] Initial test of script
	- [ ] Use this for testing Ethanol
		- [ ] Going to look at different UA and AA models (listed above)
			- [ ] **AA Models**
				- [ ] OPLS-AA
				- [ ] OPLS-AA/M 
				- [ ] GaFF/GaFF2
				- [ ] Charm
				- [ ] Amber
			- [ ] **UA Models**
				- [ ] TraPPE-UA
				- [ ] GROMOS
				- [ ] PolCA