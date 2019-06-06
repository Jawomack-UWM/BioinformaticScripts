This was the original way I ran sequnecing when I first got the position.  
Now everything is organized better into folders that allows for our shiny App "IRNAA" to perform most the analysis

Typical path for salmon -->
1) If I have already annotated the downloaded transcriptome  I will use that.  Otherwise I annotate with r package biomaRt.
biomart.R

2) I perform quantification using Looping Salmon.py.  This allows me to quantify all the files in a bash through python.
i.e. go home and let the computer work.

3. I combine all the counts from the quant.sf file outputs from salmon.py into one text file using concat_quant.py

Then these can be prepared for further downstream analysis depending on desired outcome.
If desired transcripts are converted to gene level abundance using the bioconductor package 'tximport.'
Differential expression R code can be seen in my 'RNA SEQ STAR PROJECT'

