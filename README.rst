eCGF
In silico CGF subtyping of Campylobacter using whole-genome sequence data.
Prediction accuracy of 1046 genomes: 99.682%

Getting Started
How to access and install the program:

One time:
1) Download and install Blast+:
    ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-win64.exe
2) Download and install Python:
    NOTE: This also installs pip.
    a) Go to: python.org/downloads
    b) Click download link next to python 3.6.3 (SHAN: CHECK VERSION USED)
    c) IMPORTANT: Click box: Add Python 3.6 to PATH
    d) Click Install Now
    e) Click yes
    f) Click Disable Path Length limit (doesn’t really matter but doesn’t hurt).

3) Download eCGF Program
    a) Go to my Github account releases:
        https://github.com/sfisher4/CGFPrediction/releases
    b) Click “source code zip” of the latest release (topmost release) to download the zip file of the latest release.

4) Download dependencies:
    a) Go to start menu, type powershell and open up windows powershell.
    b) Type in shell: pip install memory_profiler bioseq biopython
            NOTE: If it says pip cannot be found, then you didn't add python to your PATH when installing it,
            so reinstall and make sure you pressed the add to path button.
            If that doesn’t work: you can just replace all of your pip commands with
            :C:\Users\<username>\AppData\Local\Programs\Python\Python36-32\Scripts\pip.exe
    c) cd to wherever the CGFPrediction zip file from step 3 is located.
            i.e. Type in shell: cd ~\Downloads     (unless you moved the zip file location)
    d) Type in shell: pip install .\CGFPrediction-1.9.zip
    e) TEST: Type in shell: eCGF --help
        If this test works --> You have successfully downloaded and installed python and eCGF.
        If this test doesn't work --> Read error to see if it gives any clues on where it failed and/or review steps.

Each time you want to run eCGF:
1) Go to start menu, type powershell and open up windows powershell.
2) Type: eCGF "./path/to/folder/containing/all/fasta/genomes" "./path/to/output.csv"
NOTE: Every fasta file in the input folder will be run through the program.

** Still testing below to make sure works. Will test on Dec 6th or 7th at home **
To Add Fingerprints to the Database:
    1) Unzip CGFPrediction folder.
    2) Open cgf_pred/csvs. In this folder you will find a file called "cgf_types_fprints"
    3) Modify or replace this file ensuring formatting is not changed.
    4) Compress the CGFPrediction folder. To avoid confusion, I recommend replacing the old zip file that you
        decompressed in step 1.

Result:
A csv file containing:
    - genome names
    - eCGF binary fingerprints (2's not included, will take binary of closest match at 2's position)
    - Nearest Match	(compared to database located in file cgf_types_fprints)
    - Number of Similarities (/40)
    - Genes with Disagreement / Error Rate (avg: 0.00425)
    - Genes with 2:
    - Other cgf.type with same # of matches from 2 Case:
    - Other cgf.type with same # of matches:


**NOTES: