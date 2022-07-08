*setta il font arial*/
foreach dim in fontface fontfacesans fontfaceserif   fontfacemono {  
	graph set eps `dim' Arial
	}
global larghezzafigura "5.03"

global nero `"0 0 0 1"'
global bianco `"0 0 0 0"'
global graphregion "fcolor("${bianco}") lcolor("${bianco}") ifcolor("${bianco}") ilcolor("${bianco}")"
global plotregion "fcolor("${bianco}") lcolor("${bianco}") ifcolor("${bianco}") ilcolor("${bianco}") "
global legendstyle " region(lcolor("${bianco}") )  keygap(1) symxsize(6) size(small)"
global subtitle "subtitle(,bcolor("${bianco}") size(small))"
global ygrid "glwidth(thin) glcolor(gray*.2)"


global regionoptions "graphregion(${graphregion}`marginregion') plotregion(${plotregion} `plotregione' )"

global colorMaschi `"1 .9 0 0"'  /*blue*/
global colorFemmine `".16 1 .93 .02"' 

global color1 `".11 .49 .98 0"' /*arancio*/
global color2 `".3 1 .4 0"' /*fucsia*/
global color3 `".75 1 0 0"' /*viola*/
//global color4 `".79 .27 .08 0"' /*azzurro*/
global color5 `".53 0 .26 0"' /*celeste chiaro*/
global color6 `".24 .24 .97 0"' /*giallo*/
global color7 `".61 .19 1 .02"' /*verde */
global color4 `"0 0 0 .3"' /*grigio*/
global color8 `"0 0 0 1"' /*nero*/

global color9 `".16 1 .93 .02"' /*rosso*/
global color11 `".05 .8 1 0"' /*marrone scuro*/
global color10 `"0 .2 .5 0"' /*rosa chiaro*/
global color12 `"0 1 .7 .1"' /*magenta*/
global color13 `".06 .21 1 0"' /*ocra*/
//global color14 `".04 .85 .83 0"' /*rosso mattone*/
global color14 `"0 .35 .5 0"' /*marrone chiaro*/

global color16 `"1 0 .36 .34"' /*altro verde (giada)*/
global color17 `"0 .04 .64 .45"' /*ancora altro verde (oliva)*/


global coloraflibercept `".16 1 .93 .02"' /*rosso*/
global colorbevacizumab `".61 .19 1 .02"' /*verde */
global colorranibizumab `"1 .9 0 0"'  /*blue*/
global colordexamethasone  `".53 0 .26 0"' /*celeste chiaro*/
global colorpegaptanib `".11 .49 .98 0"' /*arancio*/

global colorgrad1 `"0 0 0 .3"' /*grigio*/
global colorgrad2 `"0 0 0 .4"' 
global colorgrad3  `"0 .2 .5 0"' /*rosa chiaro*/
global colorgrad4 `".06 .21 1 0"' /*ocra*/
global colorgrad5 `"0 .35 .5 0"' /*marrone chiaro*/
global colorgrad6 `".3 1 .4 0"' /*fucsia*/
global colorgrad7 `".75 1 0 0"' /*viola*/
global colorgrad8 `".16 1 .93 .02"' /*rosso*/
global colorgrad9 `"0 0 0 1"'
