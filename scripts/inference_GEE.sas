%macro loop(nmax);

%do modelno=0 %to &nmax;

PROC IMPORT OUT= WORK.SET
            DATAFILE= "R:\EPI\MOD\Projects\Prometheus\Work-data_scripts_programs\Graphical_Irene_new\basecase_subpopulation\basecase_subpopulation_equilibrium_nTypes5_k3_h1_pIntra0_25_c3_nObjects100000_set&modelno..csv"
            DBMS=DLM REPLACE;
     DELIMITER=',';
     GETNAMES=YES;
     DATAROW=2;
RUN;

PROC SORT DATA=SET;
  BY VAR1 RISK;
RUN;

PROC TRANSPOSE DATA=SET OUT=SWEET;
  BY VAR1 RISK;
RUN;

DATA SWEET;
  SET SWEET;
  ID = VAR1;
  TYPE = _NAME_;
  STATUS = COL1;
  KEEP ID RISK TYPE STATUS;
RUN;

ods output "Analysis Of GEE Parameter Estimates - Empirical Std Errors"=EST;
PROC GENMOD DATA=SWEET DESCENDING;
  CLASS ID TYPE STATUS;
  MODEL STATUS = RISK TYPE / LINK=LOGIT DIST=BIN;
  REPEATED SUBJECT=ID / WITHINSUBJECT=TYPE LOGOR=FULLCLUST;
RUN;
quit;
ods output close;

DATA EST_&modelno;
  SET EST;
  Model = &modelno;
  IF Parm not in ('Intercept','TYPE') THEN OUTPUT;
  KEEP Model Parm Estimate ProbZ;
RUN;

%end;
%mend;

%loop(99);

DATA EST_100;
  SET EST_0 EST_1 EST_2 EST_3 EST_4 EST_5 EST_6 EST_7 EST_8 EST_9
      EST_10 EST_11 EST_12 EST_13 EST_14 EST_15 EST_16 EST_17 EST_18 EST_19
      EST_20 EST_21 EST_22 EST_23 EST_24 EST_25 EST_26 EST_27 EST_28 EST_29
      EST_30 EST_31 EST_32 EST_33 EST_34 EST_35 EST_36 EST_37 EST_38 EST_39
      EST_40 EST_41 EST_42 EST_43 EST_44 EST_45 EST_46 EST_47 EST_48 EST_49
      EST_50 EST_51 EST_52 EST_53 EST_54 EST_55 EST_56 EST_57 EST_58 EST_59
      EST_60 EST_61 EST_62 EST_63 EST_64 EST_65 EST_66 EST_67 EST_68 EST_69
      EST_70 EST_71 EST_72 EST_73 EST_74 EST_75 EST_76 EST_77 EST_78 EST_79
      EST_80 EST_81 EST_82 EST_83 EST_84 EST_85 EST_86 EST_87 EST_88 EST_89
      EST_90 EST_91 EST_92 EST_93 EST_94 EST_95 EST_96 EST_97 EST_98 EST_99;
RUN;

PROC SORT DATA=EST_100;
  BY Model;
RUN;

PROC TRANSPOSE DATA=EST_100 OUT=TEST_100;
  BY Model;
RUN;

PROC EXPORT DATA= TEST_100
            OUTFILE= "N:\Prometheus\Manuscripts\5th paper\definitive analyses\basecase_subpopulation\combined\GEE_p estimates_N100000_c.csv"
            DBMS=CSV REPLACE;
RUN;
