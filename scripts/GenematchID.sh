#!/bin/bash
grep -w "CDS" Aspergillus_flavus.genomic.gff | cut -f9 | sed 's/Name=/\t/;s/locus_tag=/\t/' |  cut -f2,3 | sed 's/;/\t/' | cut -f1,3 | sed 's/;/\t/' | cut -f1,2 | sort | uniq > Aspergillus_flavus.GenematchID
grep -w "CDS" Aspergillus_fumigatus.genomic.gff | cut -f9 | sed 's/Name=/\t/;s/locus_tag=/\t/' |  cut -f2,3 | sed 's/;/\t/' | cut -f1,3 | sed 's/;/\t/' | cut -f1,2 | sort | uniq > Aspergillus_fumigatus.GenematchID
grep -w "CDS" Colletotrichum_gloeosporioides.genomic.gff | cut -f9 | sed 's/Name=/\t/;s/locus_tag=/\t/' |  cut -f2,3 | sed 's/;/\t/' | cut -f1,3 | sed 's/;/\t/' | cut -f1,2 | sort | uniq | grep "XP" > Colletotrichum_gloeosporioides.GenematchID
grep -w "CDS" Fusarium_oxysporum.genomic.gff | cut -f9 | sed 's/Name=/\t/;s/locus_tag=/\t/' |  cut -f2,3 | sed 's/;/\t/' | cut -f1,3 | sed 's/;/\t/' | cut -f1,2 | sort | uniq > Fusarium_oxysporum.GenematchID
grep -w "CDS" Rhizoctonia_solani.genomic.gff | cut -f9 | sed 's/Name=/\t/;s/locus_tag=/\t/' |  cut -f2,3 | sed 's/;/\t/' | cut -f1,3 | sed 's/;/\t/' | cut -f1,2 | sort | uniq > Rhizoctonia_solani.GenematchID
