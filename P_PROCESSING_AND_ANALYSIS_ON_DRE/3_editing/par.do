set more off

global dirinput "../input/ARS/g_export/Dummy tables October/"
global dirinputpoisson "./i_input_poisson/"
global dirinputpoisson "../2_analysis\pooled_input_from_DAPs/"
global dirfig "./g_figures/"
global dirpdf "./g_figures_pdf/"
global dirpng "./g_figures_png/"
global dirtab "./g_tables/"
global dircsv "./g_csv/"
global dirdtain "./g_dta_input/"
global dirdta "./g_dta/"
global dirout "./g_output/"
global dirtabtex "./g_tables_tex/"
foreach dir in dirfig dirpdf dirtab dirdtain dirdta dirout dirpng dirtabtex dircsv{
	capture mkdir ${`dir'}
	}

global datasources "ARS BIFAP_PC BIFAP_PC_HOSP PHARMO CPRD"
global datasourcesdescr "$datasources"


global nameARS "IT-ARS"
global nameBIFAP_PC "ES-BIFAP-PC"
global nameBIFAP_PC_HOSP "ES-BIFAP-PC-HOSP"
global namePHARMO "NL-PHARMO"
global nameCPRD "UK-CPRD"
	
global namesuffARS "Italy_ARS"
global namesuffBIFAP_PC "ES_BIFAP"
global namesuffBIFAP_PC_HOSP "ES_BIFAP"
global namesuffPHARMO "Netherlands-PHARMO"
global namesuffCPRD "UK_CPRD"

// global namePfizer "Pfizer"
// global nameAZ "AstraZeneca"
// global nameAstraZeneca "AstraZeneca"
// global nameModerna "Moderna"
// global nameJ&J "Janssen"
// global nameJanssen "Janssen"
// global nameUKN "Unknown"

global namePfizer "Comirnaty"
global nameAZ "Vaxzevria"
global nameAstraZeneca "Vaxzevria"
global nameModerna "Spikevax"
global nameJ&J "Janssen vaccine"
global nameJanssen "Janssen vaccine"
global nameJJ "Janssen vaccine"
global nameUKN "Unknown"

global numtab10ARS "7"
global numtab10BIFAP_PC "10"
global numtab10BIFAP_PC_HOSP "10"
global numtab10PHARMO "8"
global numtab10CPRD "9"

global months "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec"


// global events_composite "CVSTNoTP CVSTTP ArterialNoTP ArterialTP ArterialVTENoTP ArterialVTETP VTENoTP VTETP"
// global events_composite " CVSTTP  ArterialTP  ArterialVTETP  VTETP"
global events_composite "  ArterialTP VTETP ArterialVTETP "
global eventsCOVID "COVID_L1plus COVID_L2plus COVID_L3plus COVID_L4plus COVID_L5plus "
global events_narrow_broad "GBS ADEM NARCOLEPSY ACUASEARTHRITIS DM1 TP MICROANGIO HF STRCARD CAD ARR MYOCARD Myocardalone PERICARD DIC VTE TMA Hemostroke Ischstroke Sinusthrom SOCV ALI AKI GENCONV MENINGOENC  TRANSMYELITIS BP ARD ERYTH CHILBLAIN ANOSMIA ANAPHYL MISCC MIS KD DEATH SUDDENDEAT"

// global events "ACUASEARTHRITIS ADEM AKI ALI ANAPHYL ANOSMIA ARD ARR BP CAD CHILBLAIN DEATH DIC ERYTH GBS GENCONV HF Hemostroke Ischstroke KD MENINGOENC MICROANGIO MIS MYOCARD Myocardalone  PERICARD  NARCOLEPSY  SOCV STRCARD SUDDENDEAT Sinusthrom TMA TP TRANSMYELITIS VTE "

global events "GBS ADEM NARCOLEPSY ACUASEARTHRITIS DM1 TP MICROANGIO HF STRCARD CAD ARR MYOCARD Myocardalone PERICARD DIC VTE TMA Hemostroke Ischstroke Sinusthrom ${events_composite} SOCV ALI AKI GENCONV MENINGOENC  TRANSMYELITIS BP ARD ERYTH CHILBLAIN ANOSMIA ANAPHYL MISCC DEATH SUDDENDEAT ${eventsCOVID}"
global eventsnoCOVID "GBS ADEM NARCOLEPSY ACUASEARTHRITIS  TP MICROANGIO HF STRCARD CAD ARR MYOCARD Myocardalone PERICARD DIC VTE TMA Hemostroke Ischstroke Sinusthrom ${events_composite} SOCV ALI AKI GENCONV MENINGOENC  TRANSMYELITIS BP ARD ERYTH CHILBLAIN ANOSMIA ANAPHYL MISCC DEATH SUDDENDEAT"


global systems "AUTOIMMSYS CARDVASS CIRCULSYS HEPSYS NERVSYS RESPSYS SKINSYS OTHERSYS"


global eventsAUTOIMMSYS "GBS ADEM NARCOLEPSY ACUASEARTHRITIS DM1 TP"
global eventsCARDVASS "MICROANGIO HF STRCARD CAD ARR MYOCARD Myocardalone PERICARD"
global eventsCIRCULSYS "DIC VTE TMA Hemostroke Ischstroke Sinusthrom CVSTTP  ${events_composite} SOCV"
global eventsHEPSYS "ALI AKI"
global eventsNERVSYS "GENCONV MENINGOENC TRANSMYELITIS BP"
global eventsRESPSYS "ARD"
global eventsSKINSYS "ERYTH CHILBLAIN"
// global eventsOTHERSYS "ANOSMIA ANAPHYL MIS KD DEATH SUDDENDEAT ${eventsCOVID}"
global eventsOTHERSYS "ANOSMIA ANAPHYL MISCC DEATH SUDDENDEAT "

global AESIs "$events "

global AESIexcludedARS `"DM1"'
global AESIexcludedCPRD `"DM1"'
global AESIexcludedBIFAP_PC `"ANOSMIA DM1 TRANSMYELITIS"'
global AESIexcludedBIFAP_PC_HOSP `"ANOSMIA DM1 TRANSMYELITIS"'
global AESIexcludedPHARMO "ACUASEARTHRITIS ADEM AKI ALI  ARD  DIC DM1 ERYTH MICROANGIO MISCC NARCOLEPSY SOCV STRCARD  Myocardalone PERICARD  Sinusthrom TMA TRANSMYELITIS CVSTNoTP CVSTTP COVID_L4plus COVID_L5plus DEATH"

global sysnameAUTOIMMSYS "Auto-immune diseases" 
global sysnameCARDVASS "Cardiovascular system"
global sysnameCIRCULSYS "Circulatory system"
global sysnameHEPSYS "Hepato-gastrointestinal and renal system"
global sysnameNERVSYS "Nerves and central nervous system"
global sysnameRESPSYS "Respiratory system"
global sysnameSKINSYS "Skin and mucous membrane, bone and joints system"
global sysnameOTHERSYS "Other systems"

global nameACUASEARTHRITIS "Acute Aseptic Arthritis"
global nameADEM "Acute disseminated myelitis"
global nameAKI "Acute Kidney Injury"
global nameALI "Acute Liver Injury"
global nameANAPHYL "Anaphylaxis"
global nameANOSMIA "Anosmia/dysgeusia"
global nameARD "Acute respiratory distress"
global nameARR "Arrythmia"
global nameBP "Bell's palsy"
global nameCAD "Coronary artery disease"
global nameCHILBLAIN  "Chilblain"
global nameDEATH "Death"
global nameDIC "Disseminated intravascular coagulation"
global nameDM1 "Diabetes mellitus type 1"
global nameERYTH "Erythema multiforme" 
global nameGBS "Guillain Barre Syndrome"
global nameGENCONV "Generalized convulsions" 
global nameHF "Heart failure"
global nameHemostroke "Hemorragic stroke" 
global nameIschstroke "Ischemic stroke"
global nameKD "Kawasaki disease"
global nameMENINGOENC "Meningoencephalitis"
global nameMICROANGIO "Microangiopathy"
global nameMIS "Multi-inflammatory syndrome"
global nameMISCC "Multi-inflammatory syndrome"
global nameMYOCARD "Myocarditis or pericarditis"
global nameMyocardalone "Myocarditis"
global nameNARCOLEPSY "Narcolepsy"
global namePERICARD "Pericarditis"
global nameSOCV "Single Organ Cutaneous Vasculitis"
global nameSTRCARD "Stress Cardiomyopathy"
global nameSUDDENDEAT "Sudden death" 
global nameSinusthrom "Cerebral venous sinus thrombosis"
global nameTMA "Thrombotic microangiopathy"
global nameTP "Thrombocytopenia"
global nameTRANSMYELITIS "Transverse myelitis" 
global nameVTE "Venous thromboembolism"

global nameCVSTNoTP "Diagnosed or possible cerebral venous sinus thrombosis, without concurrent thrombocytopenia" 
global nameCVSTTP "Diagnosed or possible cerebral venous sinus thrombosis, with concurrent thrombocytopenia"
global nameArterialNoTP "Diagnosed or possible coronary artery disease or ischemic stroke, without concurrent thrombocytopenia" 
global nameArterialTP "Diagnosed or possible coronary artery disease or ischemic stroke, with concurrent thrombocytopenia"
global nameArterialVTENoTP "Diagnosed or possible coronary artery disease or ischemic stroke or venous thromboembolism, without concurrent thrombocytopenia" 
global nameArterialVTETP "Diagnosed or possible coronary artery disease or ischemic stroke or venous thromboembolism, with concurrent thrombocytopenia"
global nameVTENoTP "Diagnosed or possible venous thromboembolism, without concurrent thrombocytopenia"
global nameVTETP "Diagnosed or possible venous thromboembolism, with concurrent thrombocytopenia"

global nameCOVID_L1plus "Any recorded level of COVID severity"
global nameCOVID_L2plus "Level of COVID severity at least 2"
global nameCOVID_L3plus "Level of COVID severity at least 3"
global nameCOVID_L4plus "Level of COVID severity at least 4"
global nameCOVID_L5plus "COVID death"




// global birthcohorts `"<1940 1940-1949 1950-1959 1960-1969 1970-1979 1980-1989 1990+"'
global agebands `"80+ 70-79 60-69 50-59 40-49 30-39 25-29 18-24 12-17 5-11 0-4"'

global captionattrition `"Attrition diagram"'
global captioncohort_characteristics `"Cohort characteristics at start of study (1-1-2020)"'
foreach datasource of global datasourcesdescr {
	foreach dose in first second{
		global caption`dose'dose`datasource' `"Cohort characteristics at `dose' dose of COVID-19 vaccine in ${name`datasource'}"'
		}
	}
global captiondosesgaps `"COVID-19 vaccination by dose and time period between first and second dose (days)"'
global captionincident_cases `"Incident cases of AESIs in 2020 and 2021"'
foreach system of global systems{
	global captionIRvacc`system' `"${sysname`system'}. Incidence rates in vaccinated cohorts within 28 days from vaccination, and during 2020."'
	}
		
/*levels of severity of COVID mesured in each datasource*/
global levelsARS " 1 2 3 4 5 "
global levelsBIFAP_PC " 1 3 4 5 "
global levelsBIFAP_PC_HOSP " 1 3 4 5 "
global levelsCPRD " 1 2 3 4 5 "
global levelsPHARMO " 1 "
global namelevelseverity1 "Any level of severity"
global namelevelseverity2 "Level of severity 2 or worse"
global namelevelseverity3 "Level of severity 3 or worse"
global namelevelseverity4 "Level of severity 4 or worse"
global namelevelseverity5 "Death with COVID-19"

//local title = cond(`level' == 1, "Any level of severity", cond(`level' == 5,"Death with COVID-19", "Level of severity `level' or worse"))
