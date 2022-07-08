/* this step creates a pdf document incorporating tables and figures*/ 

set more off
qui do par.do

global namedoc "longpdf_220131"


/*********************************************************/
/*                                                       */
/*      DOCUMENT                                         */
/*                                                       */
/*********************************************************/
// file write `rel' `"\resizebox*{\textwidth}{!}{\input{\dirtabdoc/${anno}`nomeindic'_asl.tex}}"' _n
local namedoc "tables_and_figures"
local namedoc "$namedoc"
global dirfacsimile "./p_facsimile/"
global maketitle "\maketitle"
global tableofcontents "\tableofcontents\label{TOC}"
global listoftables "\listoftables "
global listoffigures "\listoffigures \label{TOF}"
global dirfigdoc "${dirpdf}"
global dirtabdoc "${dirtabtex}/"
// 	global dircaptdoc "../../../${dirtab}/"
global titolo_documento "Report November 2021 updated: tables and figures"
rewrite ${dirfacsimile}preambolo.tex using ${dirrep}`namedoc'.tex,replace
tempname rel
file open `rel' using ${dirrep}`namedoc'.tex, write append
file write `rel' `"\renewcommand{\arraystretch}{1.0}"' _n _n
file write `rel' `"\clearpage"' _n _n


file write `rel' `"\section{Characteristics of the cohort}"' _n _n
file write `rel' `"\subsection{Attrition diagram }"' _n _n

file write `rel' `"\begin{table}[hb]"' _n _n
file write `rel' `"\caption{${captionattrition}}"' _n _n
file write `rel' `"\label{attrition}"' _n _n
file write `rel' `"\footnotesize"' _n _n
file write `rel' `"\input{\dirtabdoc/attrition.tex}"' _n _n
file write `rel' `"\vspace{5mm}"' _n _n
file write `rel' `"\hyperlink{TOT}{\footnotesize Back to List of Tables}"' _n _n
file write `rel' `"\end{table}"' _n _n

file write `rel' `"\subsection{Characteristics of the cohort at study start}"' _n _n


capture confirm file ${dirfigdoc}PYs_overall.pdf
if _rc == 0 {
	file write `rel' `"\clearpage"' _n _n
	local caption `"Person years in the study, per week and data source"'
//		file write `rel' `"\begin{sidewaysfigure}[h]\caption{`caption'}\label{`AESI'_fig_ir_week_since_vaccination}"' _n
	file write `rel' `"\begin{figure}[h]\caption{`caption'}\label{PYs_overall}"' _n
	file write `rel' `"\includegraphics[width=.99\textwidth]{\dirfigdoc/PYs_overall.pdf}"' _n
	file write `rel' `"\vspace{5mm}"' _n _n
	file write `rel' `"\hyperlink{TOF}{\footnotesize Back to List of Figures}"' _n _n

	file write `rel' `"\end{figure}"' _n 
	file write `rel' `"\clearpage"' _n _n
	}	


file write `rel' `"\begin{table}[hb]"' _n _n
file write `rel' `"\caption{${captioncohort_characteristics}}"' _n _n
file write `rel' `"\label{cohort_characteristics}"' _n _n
file write `rel' `"\footnotesize"' _n _n
file write `rel' `"\resizebox*{\textwidth}{!}{"' _n
file write `rel' `"\input{\dirtabdoc/cohort_characteristics.tex}"' _n _n
file write `rel' `"}"' _n _n
file write `rel' `"\vspace{5mm}"' _n _n
file write `rel' `"\hyperlink{TOT}{\footnotesize Back to List of Tables}"' _n _n
file write `rel' `"\end{table}"' _n _n

file write `rel' `"\clearpage"' _n _n




file write `rel' `"\section{Vaccine exposure}"' _n _n
file write `rel' `"\subsection{Characteristics of the cohorts at first and second dose in each data source}"' _n _n

foreach datasource of global datasourcesdescr {
	foreach dose in first second{
		file write `rel' `"\begin{table}[hb]"' _n _n
		file write `rel' `"\caption{${caption`dose'dose`datasource'}}"' _n _n
		file write `rel' `"\label{`dose'dose`datasource'}"' _n _n
		file write `rel' `"\footnotesize"' _n _n
//		file write `rel' `"\resizebox*{\textwidth}{!}{"' _n
		file write `rel' `"\input{\dirtabdoc/`dose'dose`datasource'.tex}"' _n _n
//		file write `rel' `"}"' _n
		file write `rel' `"\vspace{5mm}"' _n _n
		file write `rel' `"\hyperlink{TOT}{\footnotesize Back to List of Tables}"' _n _n

		file write `rel' `"\end{table}"' _n _n
		}
	}

	

file write `rel' `"\clearpage"' _n _n

file write `rel' `"\subsection{Number of doses and gap between doses}"' _n _n

	
file write `rel' `"\begin{table}[hb]"' _n _n
file write `rel' `"\caption{${captiondosesgaps}}"' _n _n
file write `rel' `"\label{dosesgaps}"' _n _n
file write `rel' `"\footnotesize"' _n _n
file write `rel' `"\input{\dirtabdoc/dosesgaps.tex}"' _n _n
file write `rel' `"\vspace{5mm}"' _n _n
file write `rel' `"\hyperlink{TOT}{\footnotesize Back to List of Tables}"' _n _n
file write `rel' `"\end{table}"' _n _n

file write `rel' `"\clearpage"' _n _n

file write `rel' `"\subsection{Vaccine coverage per age band}"' _n _n


foreach datasource of global datasources{
	capture confirm file ${dirfigdoc}cum_coverage_`datasource'.pdf
	if _rc == 0 {
		file write `rel' `"\clearpage"' _n _n
		local caption `"Comulative weekly coverage of vaccinations in ${name`datasource'}, stratified per ageband"'
//		file write `rel' `"\begin{sidewaysfigure}[h]\caption{`caption'}\label{`AESI'_fig_ir_week_since_vaccination}"' _n
		file write `rel' `"\begin{figure}[h]\caption{`caption'}\label{cum_coverage_`datasource'}"' _n
		file write `rel' `"\includegraphics[width=.99\textwidth]{\dirfigdoc/cum_coverage_`datasource'}"' _n
		file write `rel' `"\vspace{5mm}"' _n _n
		file write `rel' `"\hyperlink{TOF}{\footnotesize Back to List of Figures}"' _n _n

		file write `rel' `"\end{figure}"' _n 
		file write `rel' `"\clearpage"' _n _n
		}	
	} 


file write `rel' `"\section{Levels of severity of COVID-19}"' _n _n
file write `rel' `"\subsection{Per calendar time in unvaccinated and vaccinated cohorts}"' _n _n
local bandlist "young old"
foreach datasource of global datasources{
	foreach band of local bandlist{
		local Age = cond("`band'" == "young","0-59","60+")
		// "${dirpdf}covidlevels_`datasource'_`band'.pdf",replace
		capture confirm file ${dirfigdoc}covidlevels_`datasource'_`band'.pdf
		if _rc == 0 {
			file write `rel' `"\clearpage"' _n _n
			local caption `"Crude weekly incidence rate  of cases of COVID-19 in ${name`datasource'}, per 100,000 person-years, by severity level, in vaccinated and unvaccinated cohorts, per calendar week since January 2020, ageband `Age'. Red crosses indicate outliers and the exact values are removed for the sake of clarity."'
	//		file write `rel' `"\begin{sidewaysfigure}[h]\caption{`caption'}\label{`AESI'_fig_ir_week_since_vaccination}"' _n
			file write `rel' `"\begin{figure}[h]\caption{`caption'}\label{covidlevels_`datasource'_`band'}"' _n
			file write `rel' `"\includegraphics[width=.99\textwidth]{\dirfigdoc/covidlevels_`datasource'_`band'.pdf}"' _n
			file write `rel' `"\vspace{5mm}"' _n _n
			file write `rel' `"\hyperlink{TOF}{\footnotesize Back to List of Figures}"' _n _n
			file write `rel' `"\end{figure}"' _n 
			file write `rel' `"\clearpage"' _n _n
			}
		}
	}

file write `rel' `"\subsection{Per time since vaccination in vaccinated cohorts}"' _n _n
foreach datasource of global datasources{
	foreach vax in Pfizer AstraZeneca Moderna Janssen{
		foreach band of local bandlist{
	
			//"${dirpdf}covidlevels_timesincevax_`vax'_`datasource'_`band'.pdf"
			local Age = cond("`band'" == "young","0-59","60+")
			capture confirm file ${dirfigdoc}covidlevels_timesincevax_`vax'_`datasource'_`band'.pdf
			if _rc == 0 {
				file write `rel' `"\clearpage"' _n _n
				local caption `"Crude weekly incidence rate of cases of COVID-19 in ${name`datasource'}, per 100,000 person-years, by levels of  severity, in cohorts aged `Age' vaccinated with ${name`vax'}, per week since vaccination."'
		//		file write `rel' `"\begin{sidewaysfigure}[h]\caption{`caption'}\label{`AESI'_fig_ir_week_since_vaccination}"' _n
				file write `rel' `"\begin{figure}[h]\caption{`caption'}\label{covidlevels_timesincevax_`vax'_`datasource'_`band'}"' _n
				file write `rel' `"\includegraphics[width=.99\textwidth]{\dirfigdoc/covidlevels_timesincevax_`vax'_`datasource'_`band'.pdf}"' _n
				file write `rel' `"\vspace{5mm}"' _n _n
				file write `rel' `"\hyperlink{TOF}{\footnotesize Back to List of Figures}"' _n _n
				file write `rel' `"\end{figure}"' _n 
				file write `rel' `"\clearpage"' _n _n
				}
			}
		}
	}
	
	
file write `rel' `"\section{Adverse events of special interest}"' _n _n
file write `rel' `"\subsection{Incident cases in 2020 and 2021}"' _n _n

file write `rel' `"\begin{table}[hb]"' _n _n
file write `rel' `"\caption{${captionincident_cases}}"' _n _n
file write `rel' `"\label{incident_cases}"' _n _n
file write `rel' `"\footnotesize"' _n _n
file write `rel' `"\input{\dirtabdoc/incident_cases.tex}"' _n _n
file write `rel' `"\vspace{5mm}"' _n _n
file write `rel' `"\hyperlink{TOT}{\footnotesize Back to List of Tables}"' _n _n
file write `rel' `"\end{table}"' _n _n

file write `rel' `"\clearpage"' _n _n



local AESIlist " $eventsnoCOVID"
foreach AESI in `AESIlist'{
	file write `rel' `"\subsection{${name`AESI'}}"' _n _n
	capture confirm file ${dirfigdoc}labelled_forest_plot_`AESI'.pdf
	if _rc == 0 {
//		file write `rel' `"\clearpage"' _n _n
		local caption "Background incidence rate in 2020 of ${name`AESI'} per 100,000 person-years with exact 95\% CI, by age band"
//		file write `rel' `"\begin{sidewaysfigure}[h]\caption{`caption'}\label{`AESI'_fig_ir_week_since_vaccination}"' _n
		file write `rel' `"\begin{figure}[h]\caption{`caption'}\label{`AESI'_fig_ir_week_since_vaccination}"' _n
		file write `rel' `"\includegraphics[width=.8\textwidth]{\dirfigdoc/labelled_forest_plot_`AESI'.pdf}"' _n
		file write `rel' `"\vspace{5mm}"' _n _n
		file write `rel' `"\hyperlink{TOF}{\footnotesize Back to List of Figuress}"' _n _n

		file write `rel' `"\end{figure}"' _n 
		file write `rel' `"\clearpage"' _n _n
		}
	local manufacturerlist "Pfizer AstraZeneca Moderna Janssen"
	local bandlist "young old"
	foreach manufacturer in `manufacturerlist'{
		foreach band in `bandlist'{
			capture confirm file "${dirfigdoc}IR_cumulative_per_week_after_vaccination_`AESI'_`manufacturer' `band'.pdf"
			if _rc == 0 {
				local Age = cond("`band'" == "young","0-59","60+")
		//		file write `rel' `"\clearpage"' _n _n
				local caption "Incidence of ${name`AESI'} per 100,000 person-years, with 95\% CI, in cohorts of age `Age' vaccinated with ${name`manufacturer'}, cumulated every week after vaccination, compared with the background incidence rate in 2020(in red)"
		//		file write `rel' `"\begin{sidewaysfigure}[h]\caption{`caption'}\label{`AESI'_fig_ir_week_since_vaccination}"' _n
				file write `rel' `"\begin{figure}[h]\caption{`caption'}\label{IR_cumulative_per_week_after_vaccination_`AESI'_`manufacturer' `band'}"' _n
				file write `rel' `"\includegraphics[width=.99\textwidth]{\dirfigdoc/IR_cumulative_per_week_after_vaccination_`AESI'_`manufacturer' `band'.pdf}"' _n
				file write `rel' `"\vspace{5mm}"' _n _n
				file write `rel' `"\hyperlink{TOF}{\footnotesize Back to List of Figuress}"' _n _n

				file write `rel' `"\end{figure}"' _n 
				file write `rel' `"\clearpage"' _n _n
				}
			}
		}

	capture confirm file ${dirfigdoc}AESI_std_`AESI'.pdf
	if _rc == 0 {
//		file write `rel' `"\clearpage"' _n _n
		local caption "Age-standardised incidence of ${name`AESI'} per 100,000 person-years, with 95\% CI, in vaccinated cohorts followed-up for 28 days, compared with the age-standardised background incidence rate (in red)"
//		file write `rel' `"\begin{sidewaysfigure}[h]\caption{`caption'}\label{`AESI'_fig_ir_week_since_vaccination}"' _n
		file write `rel' `"\begin{figure}[h]\caption{`caption'}\label{AESI_std_`AESI'}"' _n
		file write `rel' `"\includegraphics[width=.8\textwidth]{\dirfigdoc/AESI_std_`AESI'.pdf}"' _n
		file write `rel' `"\vspace{5mm}"' _n _n
		file write `rel' `"\hyperlink{TOF}{\footnotesize Back to List of Figuress}"' _n _n

		file write `rel' `"\end{figure}"' _n 
		file write `rel' `"\clearpage"' _n _n
		 }
	capture confirm file "${dirtabdoc}Table XX Crude and standardised incidence rates for `AESI'.tex"
	if _rc == 0 {
		file write `rel' `"\begin{table}[hb]"' _n _n
		file write `rel' `"\caption{Crude and standardised Background incidence rates for ${name`AESI'} per 100,000 person-years, and after vaccination followed-up for 28 days: number of expected and observed events and age-standardised incidence rate}"' _n _n
		file write `rel' `"\label{crude_and_std_IR`AESI'}"' _n _n
		file write `rel' `"\footnotesize"' _n _n
		file write `rel' `"\resizebox*{\textwidth}{!}{\input{\dirtabdoc/Table XX Crude and standardised incidence rates for `AESI'.tex}}"' _n _n
		file write `rel' `"\vspace{5mm}"' _n _n
		file write `rel' `"\hyperlink{TOT}{\footnotesize Back to List of Tables}"' _n _n
		file write `rel' `"\end{table}"' _n _n

		file write `rel' `"\clearpage"' _n _n
		}
	
	}



file close `rel'
rewrite ${dirfacsimile}piede.tex using ${dirrep}`namedoc'.tex,append



cd ${dirrep}
/* the main document is compiled via pdflatex */
! pdflatex `namedoc'.tex
! pdflatex `namedoc'.tex
/* the main document is opened via acrobat reader*/
// 	capture winexec AcroRd32 report.pdf
//cd ../
	
