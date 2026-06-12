              MACS Public Data Set Distribution Document
                       for Release P27

Date: October 18, 2019    MACS Web Site: http://statepi.jhsph.edu/macs/pdt.html

A.   General Information

   Cut-off date:  Cut-off date for this release of MACS Public Data Set
        is September 30, 2017. Data covered Visits 1 through 67 of MACS.

   Availability: Online Download

  The following files are essential to download:
        All files are ASCII text file. Open them with Wordpad or any text editor software.
        readme.txt (this file)
        codebook.zip: Contents 56 codebook files.
        data.zip: Contents 56 data files (file size of 46 MB will unzip into 781 MB).
        sasinp.zip: Contents 56 SAS input files

   Importing .dat files into STATA:  The best workaround to export the .data/input formats to Stata is a short SAS
   program.   Below is SAS code that users may adapt to save Stata datasets. Before using this code, the matching
   input and data files must be saved with identical prefix names (e.g. section4.dat and section4.inp).

        %macro exportToStata(dsn);
        data &dsn;


        Infile "\\insert file path of saved .data file here\&dsn..dat" lrecl=4000 missover;

        %include "\\ insert file path of saved .input file here\&dsn..inp";

        If _n_=1 then put version=;
        run;


        Proc export data =&dsn.  outfile="\\insert file path to save Stata dataset here\&dsn..dta" replace; run;
        %mend;

        %exporttostata(insert the file name prefix here, e.g. section4); run;


B. Data Components Included in This Release

   1. Original Cohort:

        The original 4,954 gay and bisexual men volunteered since the
        beginning of the MACS study in 1984.  They were followed up
        semi-annually. In each visit the data sets of physical examination,
        sections 2, 3, and 4 of questionnaires and laboratory results were
        generated.

   2. 1987-1990 New Recruit Cohort:

        From April 1987 through September 1991 recruits were opened to focus on
        minority and special target groups such as partners of the original
        cohort.  Version of new recruit baseline questionnaires of physical
        examination, sections 2, 3 and 4 were initially applied to this cohort.
        In their follow-up visits this cohort was tested and interviewed along
        with the original cohort.

   3. Neuropsychological Cohort:

        In 1987 MACS centers began administering neuropsychological tests
        to a subgroup of the original cohort.  The following forms were
        used to collect their testing results.  The NP wave specific
        interviews and testing was ended by Wave 10.  Beginning Visit 15
        of the main cohort, NP related forms were administered at the time
        of the participant's semi-annual visit.

        Form 7:
           Baseline questions - asked at subject's first NP Wave visit.

        PHASE 1
        Form 8:
           Neuropsychological screening battery - psychological tests.

              *criteria to move to PHASE 2 in CURRENT wave - scores of
              2 SD below the mean on 2 different tests OR a score of
              of 3 SD below the mean on 1 test.

           Neurological self-report questions - asks about headaches,
              vision problems, seizures, and strength/mobility problems,
              during past three months

                 *criteria to move to PHASE 2 in CURRENT wave - 3 or more
                 "yes" responses.

           ** meeting either the screening criteria OR self-report criteria
              will move the participant to PHASE 2 in CURRENT wave.

        Form 19: Reaction time

           Reaction time is a computer-based test which measures simple and
           complex reaction time through a variety of different tasks.  During
           waves 3 thru 10 this test was unique to the LA center.  Starting
           at visit 15 (April 1991), the other three centers began administering
           this test to their NP participants.  Results of the reaction time
           testing from all centers is compiled in LA.

        PHASE 2
        Form 9:
           Neuropsychological tests - further psychological testing.

              *criteria to undergo PHASE 2 again in NEXT wave - scores of
              1.5 SD below the mean on 2 different tests OR a score of
              2.5 SD below the mean on 1 test.

        Form 10:
           Neurological exam - performed by a neurologist.

              *criteria to undergo PH 2 again in NEXT wave - an answer of
              2 or 3 to variable NEURO_K(#) in Form 10.


           ** meeting either the Form 9 criteria OR Form 10 criteria
              will cause the participant to undergo PHASE 2 in the NEXT wave.


        DIAGNOSTIC TESTS
        Participants are asked to undergo these tests at random or based on the
        results of their PH 2 tests.

           1). NCV/EMG - electrodiagnostic tests

           2). MRI - magnetic resonance imaging

           3). LP - lumbar puncture to collect cerebrospinal fluid

           4). EEG - electro-encephalograph (ceased this test in 7/88)

           * Criteria to undergo PHASE #2 in NEXT wave - if any of above tests are
           judged by a neurologist to be abnormal.

   4. 2001-2003 New Recruit Cohort:

        A third enrollment of 1,351 men took place
        between October 2001 and August 2003.  This third cohort augments
        research efforts in the long term benefits and adverse effects of
        therapy. Started visit 36.5.


   5. Health Service Utilization and Economics:

     In 1990, questions related health service utilization and economics
     related to HIV-1 disease were developed and piloted during visit 13.
     These questions were finalized and incorporated in the regular visit
     questionnaire starting with MACS visit 14.


   6. Participants With Aids (PWA)

     PWA (Participants With AIDS) Form was developed to be
     used to closely follow PWAs and those participants with
     less than 200 T-cells. A phone call was attempted
     at least every three months to PWAs. If the
     participant does not complete a full Section 4 during the
     regulat visit, this short form will be used as
     an alternative where we can collect hospitalization, medication,
     and insurance information in a very brief format.

   7. Quality of Life (QOL)

     36-Item Health Survey: quality of life (introduced at Visit 21)

   8. Drug Form1: anti-virals (introduced at Visit 13)

   9. Drug Fomr2: non-antivirals (introduced at Visit 13)

  10. Antiviral Medicaion Adherence:  Antiviral medication adherence
      (introduced at Visit 30)

  11. Clinical outcome: Included reported events of AIDS, cancer, diagnosis or death.

  12. Men's Attitude Survey: Started  Visit 31.  Administered intermittingly.

  13. Timed Walking and Hand Grip Assessments: Started Visit 43.

  14. 2010 recruits: Seronegative (controls), incident sero-converters, ART-nagatives
      and experienced.

  15. Medical Abstracts: Retrospective Medical Record Abstraction

C. Associated Questionnaire Forms and Documentation

   Forms: (Scanned and converted PDF format on MACS website)

   All MACS questionnaire forms are available from:

   http://statepi.jhsph.edu/macs/forms.html

   Other documentation:

   MACS directory of investigators and archives of publication
   can be found in MACS web site (http://statepi.jhsph.edu/macs).

Your ideas on research of the MACS data and presentation of
of your research findings are welcomed to share with the MACS
investigators.  Please contact:

             Dr. Gypsyamber D'Souza (Amber)PHD
             Principal Investigator, CAMACS
             Department of Epidemiology
             Bloomberg School of Public Health
             Johns Hopkins University
             615 North Wolfe Street
             Baltimore, MD 21205

             Phone: (410) 955-4320
             FAX:   (410) 955-7587
             Email: gdsouza2@jhu.edu

D.  File Content of the Released CD-ROM

Dir  File Name               File Description
---- ------------     ---------------------------------------------
\
     readme.txt    -- This readme file.

\codebook

     cancer.cdb    -- Codebook file for cancer section of outcome form
     drugadh.cdb   -- Codebook file for Antiviral Medicaion Adherence
     drugf1.cdb    -- Codebook file for Drug Form 1
     drugf2.cdb    -- Codebook file for Drug Form 2
     frail.cdb     -- Codebook file for Timed Walking and Hand Grip Assessments
     hivstats.cdb  -- Codebook file for HIV status
     idal.dat      -- Codebook file for IADL
     lab_rslt.cdb  -- Codebook file for Laboratory Results
     macsid.cdb    -- Codebook file for Participant IDs with Recruit status and DOB year
     mas.cdb       -- Codebook file for Men's Attitute Survey
     med_abst.cdb  -- Codebook file for Retrospetive medical record abstraction
     microglb.cdb  -- Codebook file for B2-Microglobulin Test
     neoptrin.cdb  -- Codebook file for Neopterin Test
     npf08.cdb     -- Codebook file for NP Form 08, Began Visit 15
     npf10.cdb     -- Codebook file for NP Form 10, Began Visit 15
     npf18.cdb     -- Codebook file for NP Form 18, Began Visit 16
     npfrt.cdb     -- Codebook file for NP Form 19, Reaction time
     outcome.cdb   -- Codebook file for events of clinical system and endpoint
     phy_exam.cdb  -- Codebook file for Physical Examination
     pwa.cdb       -- Codebook file for Participant with AIDS
     qol.cdb       -- Codebook file for Quality of Life
     section2.cdb  -- Codebook file for Sections 2 and 3
     section4.cdb  -- Codebook file for Section 4
     v13herc.cdb   -- Codebook file for Health Utilization and Economic, Visit 13
     w00f07.cdb    -- Codebook file for Neuropsychology, Form 07, Baseline
     w00f12.cdb    -- Codebook file for Neuropsychology, Form 12, Baseline
     w01f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 01
     w01f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 01
     w01f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 01
     w02f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 02
     w02f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 02
     w02f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 02
     w03f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 03
     w03f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 03
     w03f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 03
     w04f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 04
     w04f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 04
     w04f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 04
     w05f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 05
     w05f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 05
     w05f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 05
     w06f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 06
     w06f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 06
     w06f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 06
     w07f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 07
     w07f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 07
     w07f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 07
     w08f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 08
     w08f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 08
     w08f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 08
     w09f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 09
     w09f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 09
     w09f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 09
     w10f08.cdb    -- Codebook file for Neuropsychology, Form 08, Wave 10
     w10f09.cdb    -- Codebook file for Neuropsychology, Form 09, Wave 10
     w10f10.cdb    -- Codebook file for Neuropsychology, Form 10, Wave 10
                                                                               No. of
\data                                                                          Record
                                                                              -------
     cancer.dat    -- Data file for cancer section of outcome form               1272
     drugadh.dat   -- Data file for Antiviral Medicaion Adherence               30314
     drugf1.dat    -- Data file for Drug Form 1                                102807
     drugf2.dat    -- Data file for Drug Form 2                                 23309
     frail.dat     -- Data file for Timed Walking and Hand Grip Assessments     39665
     hivstats.dat  -- Data file for HIV status                                   7343
     idal.dat      -- Data file for IADL                                        29082
     lab_rslt.dat  -- Data file for Laboratory Results                         143263
     macsid.dat    -- Data file for Participant IDs with Recruit status & DOB    7343
     mas.dat       -- Data file for Men's Attitute Survey                       14510
     med_abst.dat  -- Data file for Retrospetive medical record abstraction       874
     microglb.dat  -- Data file for B2-Microglobulin Test                       12404
     neoptrin.dat  -- Data file for Neopterin Test                              12320
     npf08.dat     -- Data file for NP Form 08, Began Visit 15                  24444
     npf10.dat     -- Data file for NP Form 10, Began Visit 15                   2334
     npf18.dat     -- Data file for NP Form 18, Began Visit 16                  72112
     npfrt.dat     -- Data file for NP Form 19, Reaction time                   22982
     outcome.dat   -- Data file for events of clinical system and endpoint       3752
     phy_exam.dat  -- Data file for Physical Examination                       138829
     pwa.dat       -- Data file for Participant with AIDS                        2352
     qol.dat       -- Data file for Quality of Life                             72805
     section2.dat  -- Data file for Sections 2 and 3                           145303
     section4.dat  -- Data file for Section 4                                  152443
     v13herc.dat   -- Data file for Health Utilization and Economic, Visit 13    1764
     w00f07.dat    -- Data file for Neuropsychology, Form 07, Baseline           4007
     w00f12.dat    -- Data file for Neuropsychology, Form 12, Baseline            271
     w01f08.dat    -- Data file for Neuropsychology, Form 08, Wave 01             566
     w01f09.dat    -- Data file for Neuropsychology, Form 09, Wave 01             153
     w01f10.dat    -- Data file for Neuropsychology, Form 10, Wave 01             149
     w02f08.dat    -- Data file for Neuropsychology, Form 08, Wave 02             325
     w02f09.dat    -- Data file for Neuropsychology, Form 09, Wave 02              95
     w02f10.dat    -- Data file for Neuropsychology, Form 10, Wave 02              92
     w03f08.dat    -- Data file for Neuropsychology, Form 08, Wave 03            1214
     w03f09.dat    -- Data file for Neuropsychology, Form 09, Wave 03             144
     w03f10.dat    -- Data file for Neuropsychology, Form 10, Wave 03              80
     w04f08.dat    -- Data file for Neuropsychology, Form 08, Wave 04            1361
     w04f09.dat    -- Data file for Neuropsychology, Form 09, Wave 04             232
     w04f10.dat    -- Data file for Neuropsychology, Form 10, Wave 04             197
     w05f08.dat    -- Data file for Neuropsychology, Form 08, Wave 05            1160
     w05f09.dat    -- Data file for Neuropsychology, Form 09, Wave 05             199
     w05f10.dat    -- Data file for Neuropsychology, Form 10, Wave 05             162
     w06f08.dat    -- Data file for Neuropsychology, Form 08, Wave 06            1187
     w06f09.dat    -- Data file for Neuropsychology, Form 09, Wave 06             253
     w06f10.dat    -- Data file for Neuropsychology, Form 10, Wave 06             261
     w07f08.dat    -- Data file for Neuropsychology, Form 08, Wave 07            1097
     w07f09.dat    -- Data file for Neuropsychology, Form 09, Wave 07             318
     w07f10.dat    -- Data file for Neuropsychology, Form 10, Wave 07             315
     w08f08.dat    -- Data file for Neuropsychology, Form 08, Wave 08            1144
     w08f09.dat    -- Data file for Neuropsychology, Form 09, Wave 08             360
     w08f10.dat    -- Data file for Neuropsychology, Form 10, Wave 08             325
     w09f08.dat    -- Data file for Neuropsychology, Form 08, Wave 09            1570
     w09f09.dat    -- Data file for Neuropsychology, Form 09, Wave 09             349
     w09f10.dat    -- Data file for Neuropsychology, Form 10, Wave 09             326
     w10f08.dat    -- Data file for Neuropsychology, Form 08, Wave 10            1155
     w10f09.dat    -- Data file for Neuropsychology, Form 09, Wave 10             274
     w10f10.dat    -- Data file for Neuropsychology, Form 10, Wave 10             260

\form

     See Website of http://statepi.jhsph.edu/macs/forms.html for the original
     questionnaire forms that were scanned and coverted into PDF format.

\sasinp

     cancer.inp    -- SAS input for cancer section of outcome form
     drugadh.inp   -- SAS input for Antiviral Medicaion Adherence
     drugf1.inp    -- SAS input for Drug Form 1
     drugf2.inp    -- SAS input for Drug Form 2
     frail.inp     -- SAS input file for Timed Walking and Hand Grip Assessments
     hivstats.inp  -- SAS input for HIV status
     idal.inp      -- SAS input for IADL
     lab_rslt.inp  -- SAS input for Laboratory Results
     macsid.inp    -- SAS input for Participant IDs with Recruit status and DOB Year
     mas.inp       -- SAS input file for Men's Attitute Survey
     med_abst.inp  -- SAS input file for Retrospetive medical record abstraction
     microglb.inp  -- SAS input for B2-Microglobulin Test
     neoptrin.inp  -- SAS input for Neopterin Test
     npf08.inp     -- SAS input for NP Form 08, Began Visit 15
     npf10.inp     -- SAS input for NP Form 10, Began Visit 15
     npf18.inp     -- SAS input for NP Form 18, Began Visit 16
     npfrt.inp     -- SAS input for NP Form 19, Reaction time
     outcome.inp   -- SAS input for events of clinical system and endpoint
     phy_exam.inp  -- SAS input for Physical Examination
     pwa.inp       -- SAS input for Participant with AIDS
     qol.cdb       -- SAS input for Quality of Life
     section2.inp  -- SAS input for Sections 2 and 3
     section4.inp  -- SAS input for Section 4
     v13herc.inp   -- SAS input for Health Utilization and Economic, Visit 13
     w00f07.inp    -- SAS input for Neuropsychology, Form 07, Baseline
     w00f12.inp    -- SAS input for Neuropsychology, Form 12, Baseline
     w01f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 01
     w01f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 01
     w01f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 01
     w02f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 02
     w02f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 02
     w02f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 02
     w03f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 03
     w03f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 03
     w03f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 03
     w04f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 04
     w04f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 04
     w04f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 04
     w05f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 05
     w05f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 05
     w05f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 05
     w06f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 06
     w06f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 06
     w06f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 06
     w07f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 07
     w07f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 07
     w07f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 07
     w08f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 08
     w08f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 08
     w08f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 08
     w09f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 09
     w09f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 09
     w09f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 09
     w10f08.inp    -- SAS input for Neuropsychology, Form 08, Wave 10
     w10f09.inp    -- SAS input for Neuropsychology, Form 09, Wave 10
     w10f10.inp    -- SAS input for Neuropsychology, Form 10, Wave 10
