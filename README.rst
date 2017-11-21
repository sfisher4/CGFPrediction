eCGF
In silico CGF subtyping of Campylobacter using whole-genome sequence data.
Prediction accuracy of 1046 genomes: 99.682%

Getting Started
How to access and install the program:

1) Download the latest release on my github page:
        https://github.com/sfisher4/CGFPrediction/releases
2) Install
    pip install --upgrade ~/location_of_zip_file_downloaded_in_step_1

How to run:
eCGF

Getting help:
eCGF --help

To Add Fingerprints to the Database:
    1) Unzip CGFPrediction folder.
    2) Open cgf_pred/csvs. In this folder you will find a file called "cgf_types_fprints"
    3) Modify or replace this file ensuring formatting is not changed.
    4) Compress the CGFPrediction folder. To avoid confusion, I recommend replacing the old zip file that you
        decompressed in step 1.

Result
 - binary gene presence/absence data

**NOTES:

If you receive more fingerprints in the database, this is how you can add it to the program:
 --> must be the same order of genes that is in the program?????