#!/bin/bash

grep ACAT MASTER_indep_aggs | awk '$2<8.71e-9' | cut -f1 | sed 's/\..*//' | sed 's/\_ENST.*//' | sort | uniq > ACAT_merged_hits_indep
grep BURDEN MASTER_indep_aggs | awk '$2<8.71e-9' | cut -f1 | sed 's/\..*//' | sed 's/\_ENST.*//' | sort | uniq > BURDEN_merged_hits_indep

awk '$2 < 8.71e-8' ACAT_ordered_tophits | grep -v singleton > ACAT_sig_hits
awk '$2 < 8.71e-9' BURDEN_ordered_tophits | grep -v singleton > BURDEN_sig_hits

cut -f1 ACAT_sig_hits | sed 's/\..*//' | sed 's/\_ENST.*//' | sort | uniq > ACAT_merged_hits
cut -f1 BURDEN_sig_hits | sed 's/\..*//' | sed 's/\_ENST.*//' | sort | uniq > BURDEN_merged_hits

cut -f1 -d'_' BURDEN_merged_hits_indep | sort | uniq -c
cut -f1 -d'_' BURDEN_merged_hits | sort | uniq -c
cut -f1 -d'_' ACAT_merged_hits_indep | sort | uniq -c
cut -f1 -d'_' ACAT_merged_hits | sort | uniq -c
