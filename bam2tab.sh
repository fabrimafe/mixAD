#!/bin/bash

samtools mpileup $1 | tr '[:lower:]' '[:upper:]' | sed 's/\$//g' | grep -v '[+-]' | sed 's/\^.//g' | 
awk -v OFS='\t' '{nA=0;nC=0;nT=0;nG=0; for (i=1;i<=length($5);i++){
if (substr($5,i,1)=="A"){nA=nA+1};
if (substr($5,i,1)=="C"){nC=nC+1};
if (substr($5,i,1)=="G"){nG=nG+1};
if (substr($5,i,1)=="T"){nT=nT+1};
};
if ( ! ((nC!=0 && nT!=0 ) || (nA!=0 && nG!=0 ) )) {print $2,nA,nG,nC,nT};
}' 
