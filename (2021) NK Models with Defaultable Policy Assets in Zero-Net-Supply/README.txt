% Amaral(2021): New-Keynesian Models with Defaultable Policy Assets in Zero-Net-Supply

% What does one need to replicate the paper?

a) MatLab (tested on 2020b Update 5)
b) Dynare (tested on 4.6.4)

% How does one can run the code?

1) Copy all files and subfolders in the Paper's folder in GitHub to a local folder in your computer
2) In MatLab, open the file mainScript.m
3) In the first section of the code, change pathImages to an existing path address in your computer
4) In the section "Run Dynare", set the variable flexiblePrices = false
5) In MatLab, open the file Closed Economy\risky_ClosedEconomy.mod (this is the Dynare file), and set @#define flexiblePrices = false
6) Run the section "Run dynare" and then the section "Collect impulse response functions"
7) Redo steps 4, 5 and 6, in this order, but now setting flexiblePrices = true and @#define flexiblePrices = true
8) Run the remaining sections of the code in order


