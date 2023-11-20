#!/bin/bash
#set -eu

build_dbsnp () {
#   Purpose: Retrieve the dbSNP data and store the relevant columns in the database file "snp156.db" (156 for version 156)
#   Output: *A file called snp156.db*

    wget -b https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.40.gz
    echo "build_dbsnp: Downloaded GCF_000001405.40.gz from hgdownload.soe.ucsc.edu/goldenPath/hg38/database/." > /dev/stderr
    awk '{ print $3, $1, $2, $5 }' <(zcat GCF_000001405.40.gz | tail -n +39) > snp156.tmp #select relevant columns from data
    echo "build_dbsnp: Extracted relevant information to snp156.tmp." > /dev/stderr
    sqlite3 snp156.db "CREATE TABLE snp156(rsid text NOT NULL, chrom text NOT NULL, pos INTEGER NOT NULL, alleles text NOT NULL);"
    sqlite3 -separator ' ' snp156.db '.import "snp156.tmp" "snp156"' #make the database with the selected columns
    sqlite3 -separator ' ' snp156.db 'CREATE INDEX idx_rsid ON snp156(rsid);' #the columns are: rs identifier, chromomosme, position, allele possibilities
    chmod 444 snp156.db
    rm snp156.tmp
    echo "build_dbsnp: dbSNP database written to snp156.db." > /dev/stderr
}


pop_selection () {


inputs=()


while IFS= read -r input; do
    # Check if the input is empty or if the user pressed Ctrl+D
    if [[ -z "$input" ]]; then
        break
    fi

    # Verify if the input is present in the merged_file.txt
    if ! grep -q -x "$input" "/scratch/richards/ethan.kreuzer/data/ancestry/valid_ancestry_codes.txt"; then
        echo "Invalid input: '$input' does not exist in the ancestry file."
        exit 1  # Exit the script with an error code
    fi

    # Add the valid input to the array
    inputs+=("$input")
done < ${ANCESTRY_SELECTION}

# Print the collected inputs


tmp_file="${SNAPTMP}/ancestry.fam.tmp"
output_file="${SNAPTMP}/ancestry.fam"
for ancestry in "${inputs[@]}"; do

    file="/scratch/richards/ethan.kreuzer/data/ancestry/${ancestry}_samples.txt~"

    cat "$file" >> "$tmp_file"

    # Add your desired operations for each element here
done

# Keep only the unique lines in the output file
sort -u "$tmp_file" > "$output_file"

# Remove the tmp_file
rm "$tmp_file"

}



init () {
#    Purpose: Helper function. It saves the SNPs in the outcome_gwas into a file and outputs the query SNPs to be processed by match()
#    Parameters: 
#             
#             1. A list of Query_SNPs (results of an exposure GWAS)
#             2. A list of Outcome_SNPs (results of an outcome GWAS)
#
#    Output: 
#    
#    rs2476601 . . . .
#    rs2596546 . . . .
#    rs2603139 . . . .
#    rs2723834 . . . .

    #cat <(zcat -f -- "${GWAS_F}") > ${SNAPTMP}/SNAP.GWAS.init #store the given outcome_snps in a file called SNAP.GWAS.init
    cat <(zcat -f -- "${GWAS_F}") | awk 'BEGIN {OFS="\n"} {if (!/^rs/) {sub(/^(chr)?/, "chr"); gsub("_", ":");} print}' > "${SNAPTMP}/SNAP.GWAS.init"

    sort "${QSNP_F}" | uniq -d > ${OUTPUT_FILE}.duplicates.txt #sort query snps alphabetically and write duplicate lines to duplicates.txt
    sort "${QSNP_F}" | uniq | awk 'BEGIN {OFS="\n"} {if (!/^rs/) {sub(/^(chr)?/, "chr"); gsub("_", ":");} print}' | awk '{ print $1,".",".",".","." }' - #Query SNP ID's followed by 4 lines. 
    echo "init: $(cat ${OUTPUT_FILE}.duplicates.txt | wc -l) duplicated SNP identifiers removed, written to duplicates.txt." > /dev/stderr #Reports quantity of duplicate SNP's
}

match () { 
#   Purpose: Identify which query SNPs are explicitely in the outcome_snps
#   Output:

#   rs2476601 rs2476601 1 0 match
#   rs2596546 rs2596546 1 0 match
#   rs2603139 rs2603139 1 0 match
#   rs2723834 . . . .

    #awk 'BEGIN {OFS="\n"} {if (!/^rs/) {sub(/^(chr)?/, "chr"); gsub("_", ":");} print}' </dev/stdin > "${SNAPTMP}/SNAP.TMP.match"

    #awk 'BEGIN {OFS="\n"} {if (!/^rs/) {sub(/^(chr)?/, "chr"); gsub("_", ":");} print}' ${SNAPTMP}/SNAP.GWAS.init > ${SNAPTMP}/SNAP.GWAS.init.tmp

    #chmod 666 ${SNAPTMP}/SNAP.GWAS.init.tmp ${SNAPTMP}/SNAP.TMP.match

    #awk 'FNR==NR { m[$1]=1; next } $2=="." && ($1 in m) { $2=$1;$3=1;$4=0;$5="match";c+=1 } { print } END { print "match: Matched",c,"SNP by exact match." > "/dev/stderr" }' \
    
    #${SNAPTMP}/SNAP.GWAS.init.tmp \
    
    #${SNAPTMP}/SNAP.TMP.match
    cat </dev/stdin > ${SNAPTMP}/SNAP.TMP.match #read output from match() or update()  into SNAP.TMP.match 

     #awk 'FNR==NR { m[$1]=1; next } $2=="." && ($1 in m) { $2=$1;$3=1;$4=0;$5="match";c+=1; print } !($2=="." && ($1 in m)) { $2=".";$3="."; $4="."; $5="."; print } END { print "match: Matched",c,"SNP by exact match." > "/dev/stderr" }' \
      #  ${SNAPTMP}/SNAP.GWAS.init \
       # ${SNAPTMP}/SNAP.TMP.match

  awk 'NR == FNR { gwas_init[$1] = 1; next }
     {
         if ($2 == ".") {
             if ($1 in gwas_init) {
                 print $1, $1, ".", ".", "match"; c += 1;
             } else {
                 print $1, ".", ".", ".", ".";
             }
         } else {
             print $1, $2, $3, $4, $5;
         }
     }
     END {
         print "match: Matched", c, "SNP by exact match." > "/dev/stderr"
     }' ${SNAPTMP}/SNAP.GWAS.init ${SNAPTMP}/SNAP.TMP.match


}

position () {
#   Purpose: It is possbile the query SNPs are RSid's but the outcome SNPs are identified by chromsome position (such as "12:12055413_TATTA_T" for example).
#            The match() function will fail to report a match in such cases, even though in actuality they refer to the same SNP. This functions looks up
#            the chromsome location of the unmatched SNP RS_id's to attempt to match on the outcome SNPs by chromosme location.
#   Output:
#           
#   rs2603139 rs2603139 1 0 match
#   rs2723834 12:12055413_TATTA_T 1 0 position
#   rs2744944 rs2744944 1 0 match
#   rs2745803 rs2745803 1 0 match

    cat </dev/stdin > ${SNAPTMP}/SNAP.TMP.position
    echo "position: Started search by genomic position, please wait." > /dev/stderr
    awk '$1 !~ /rs/ { match($1, "(chr)?([0-9]+)[^0-9]+([0-9]+)", ary); print $1, "chr"ary[2]":"ary[3]":"$2 }' ${SNAPTMP}/SNAP.GWAS.init | sed 's/:$//' > ${SNAPTMP}/SNAP.GWAS.position
    if [[ -s ${SNAPTMP}/SNAP.GWAS.position ]]; then
        ## Search dbSNP by rsid, return genomic position, etc.
        while read snp; do
  output=$(sqlite3 -separator ' ' "${UCSC_DBSNP}" "select * from snp156 where rsid=\"$snp\" LIMIT 1")
  chromosome=$(awk '{print $2}' <<< "$output")
  
  # Remove everything after the first dot
  chromosome=${chromosome%%.*}
  
  # Remove everything before the underscore
  chromosome=${chromosome#*_}
  
  # Remove leading zeros
  chromosome=$(echo "$chromosome" | sed 's/^0*//')
  
  # Modify the second field ($2) with the updated chromosome value
  modified_output=$(awk -v OFS=' ' -v chr="$chromosome" '{$2=chr} 1' <<< "$output")
  
  echo "$modified_output"
done < <(awk '$2=="." { print $1 }' ${SNAPTMP}/SNAP.TMP.position) | awk '{ chr=$2; gsub("chr", "", chr); print $1,$2,$3,$4,"chr"chr":"$3 }' > ${SNAPTMP}/SNAP.DBSNP.position
        
	awk 'FNR==NR { m[$2]=$1; next } $5 in m { print $0, m[$5] }' ${SNAPTMP}/SNAP.GWAS.position ${SNAPTMP}/SNAP.DBSNP.position > ${SNAPTMP}/SNAP.DBSNP.GWAS.position

	if [[ -s ${SNAPTMP}/SNAP.DBSNP.GWAS.position ]]; then
            awk 'BEGIN { c=0 } FNR==NR { m[$1]=$6; next } $1 in m { $2=m[$1];$3=1;$4=0;$5="position";c+=1} { print } END { print "position: Matched",c,"SNPs by position." > "/dev/stderr"}' ${SNAPTMP}/SNAP.DBSNP.GWAS.position ${SNAPTMP}/SNAP.TMP.position -
        else
            echo "position: Warning -- Analysis was not performed. The target file does not contain SNP identifiers with positional information." > /dev/stderr
            cat ${SNAPTMP}/SNAP.TMP.position
        fi
    else
        echo "position: Warning -- Analysis was not performed. The GWAS file does not contain SNP identifiers with positional information." > /dev/stderr
        cat ${SNAPTMP}/SNAP.TMP.position
    fi

}

proxy () {

#   Purpose: This is the snappy output. It returns proxys for the unmatched query SNPs by finding outcome SNPs that are in LD with the unmatched query SNP.
#   Parameters:             
#	    1. Path to a directory that contains a Reference Genome to search for proxys
#           2. Path to a file that contains the ID's of the Ethnic Population to use when finding proxys
#           3. An "r2" threshold value.
#           4. The number of proxys to return
#
#   Output:
#           rs2248372 rs2248372 1 0 match
#           rs2249059 rs78802957 1 46563 rs114028162 1 46928 rs142580331 1 49109 rs113520162 1 49209 proxy
#           rs229540 rs229540 1 0 match
#           rs2402240 rs2402240 1 0 match
#           rs244686 rs244686 1 0 match
#           rs2476601 rs2476601 1 0 match
#           rs2596546 rs2596546 1 0 match
#           rs2603139 rs2603139 1 0 match
#           rs2723834 12:12055413_TATTA_T 1 0 position
    
    cat </dev/stdin > ${SNAPTMP}/SNAP.TMP.proxy #read output from position() into proxy
    
    awk '$2=="." { print $1 }' ${SNAPTMP}/SNAP.TMP.proxy > ${SNAPTMP}/SNAP.input.proxy.tmp
    
    sort ${SNAPTMP}/SNAP.input.proxy.tmp | uniq | awk -F ':' '{if ($1 !~ /^rs/) {print $1":"$2} else {print $1}}' > ${SNAPTMP}/SNAP.input.proxy #write all the unmatched SNPs into a file called SNAP.input.proxy

    set +e #continue executing commands even if some fail (this is expected behavior)
    



#CODE TO DECREASE THE FOR LOOP

# Query DBSNP to get chromosomes you need to visit


while IFS= read -r snp; do

    if [[ "${snp}" = rs* ]]; then

        chromosome=$(sqlite3 /project/richards/ethan.kreuzer/snp156.db "SELECT chrom FROM snp156 WHERE rsid = '$snp' LIMIT 1;")

        chromosome=${chromosome%%.*}  # Remove everything after the first dot

        chromosome=${chromosome#*_}  # Remove everything before the underscore

        chromosome=$(echo "$chromosome" | sed 's/^0*//')  # Remove leading zeros

        if [[ -n $chromosome ]]; then

            if [[ ${#chromosome} -ge 2 ]]; then 

	       formatted_chromosome=$(printf "%02d" "$chromosome")
        
	    else
           
	       formatted_chromosome=$(printf "%d" "$chromosome")
            fi

            chromosomes+=("$formatted_chromosome")
        fi

    else
        
	chromosome=${snp#chr}

	chromosome=${chromosome%%:*}

	chromosomes+=("$chromosome")

	fi
    
   done < ${SNAPTMP}/SNAP.input.proxy


#CODE TO DECREASE THE LOOP ITERATIONS

sorted_chromosomes=($(printf "%s\n" "${chromosomes[@]}" | sort -u))


# Define a function to execute the PLINK command
execute_plink() {
    local chrom="$1"

    echo "proxy: Running PLINK LD analysis for chromosome ${chrom}." > /dev/stderr

    local out_prefix="${SNAPTMP}/SNAP.${chrom}.proxy"

    /scratch/richards/yiheng.chen/Plink1.9/plink --bfile "${PLINK_REF_PANEL}/${chrom}.final_id_len18" \
        --r2 --ld-window ${PLINK_WINDOW} --ld-window-r2 ${PLINK_MIN_R2} \
        --keep "${SNAPTMP}/ancestry.fam" \
        --out "${out_prefix}" \
        --ld-snp-list "${SNAPTMP}/SNAP.input.proxy" \
	--memory 10000 > /dev/null
    retVal=$?
    if [ $retVal -eq 0 ]; then
        perl -pi -e "s/[ \t]+/ /g;" "${out_prefix}.ld"
        perl -pi -e "s/^[ \t]+//g;" "${out_prefix}.ld"
        perl -pi -e "s/[ \t]+$//g;" "${out_prefix}.ld"
    else
        echo "Error: PLINK analysis for chromosome ${chrom} failed." > /dev/stderr
    fi
}


# Iterate over the sorted chromosomes and run execute_plink function in the background
for chrom in "${sorted_chromosomes[@]}"; do
    execute_plink "$chrom" 
done

# Wait for all background processes to finish
wait



#set -e

## Iterate over the input files
for file in ${SNAPTMP}/SNAP.*.proxy.ld; do
    output_file="${file}.gwas" # Create the corresponding output file name

    # Run the command for the current input file and redirect the output to the output file
    awk 'FNR==NR { m[$1]=1; next } $6 in m { print }' ${SNAPTMP}/SNAP.GWAS.init <(tail -q -n +2 "$file") > "$output_file"
done

    #Filter each proxy result file to contain only SNPs already present in the outcome SNPs. Combine each filtered file into a single file.  
    

    for file in ${SNAPTMP}/SNAP.*.proxy.ld.gwas; do
      
      output_file="${file}.bestproxy"
      
      awk '{ d=$5-$2; if (d < 0) { d=$2-$5 }; print $0, d }' "$file" | sort -k 3,3 -k7,7nr -k8,8n > "$output_file"
  
   done

   




snp_db_file=/project/richards/ethan.kreuzer/snp156.db

# Function to check if an rsid is biallelic in the SNP db
function is_biallelic() {

  local rsid=$1

  if [[ $rsid != rs* ]]; then
    echo 1  # Not an rsID, return 1
    return
  fi

  local result=$(sqlite3 "$UCSC_DBSNP" "SELECT * FROM snp156 WHERE alleles NOT LIKE '%,%' AND rsid = '$rsid';")
  if [[ -n $result ]]; then
    echo 1  # Biallelic SNP
  else
    echo 0  # Not a biallelic SNP
  fi
}

# Loop over the files


for file in ${SNAPTMP}/SNAP.*.proxy.ld.gwas.bestproxy; do
  output_file="${file}.out"

  current_str=""
  counter=0

  while read -r line; do
    field3=$(echo "$line" | awk '{print $3}')
    field6=$(echo "$line" | awk '{print $6}')
    if [[ "$field3" != "$current_str" ]]; then
        current_str="$field3"
        counter=0
    fi

    if [[ "$counter" -lt ${NUM_PROXYS} ]]; then
        echo "$line" >> "$output_file"
        ((counter++))
    fi

  done < "$file"
done


   #for file in ${SNAPTMP}/SNAP.*.proxy.ld.gwas.bestproxy.out; do
  #cat "$file" >> "${SNAPTMP}/SNAP.proxy.final"
#done
  
  for file in ${SNAPTMP}/SNAP.*.proxy.ld.gwas.bestproxy.out; do
  while IFS= read -r line; do
    # Get the 6th field (assuming fields are space-separated, adjust the delimiter accordingly)
    field6=$(echo "$line" | awk '{print $6}')

    # Determine if the 6th field is biallelic using the is_biallelic function
    is_biallelic_result=$(is_biallelic "$field6")

    # Modify the line based on the is_biallelic result
    if [[ $is_biallelic_result -eq 1 ]]; then
      # If biallelic, echo the entire line as is to the final file
      echo "$line" >> "${SNAPTMP}/SNAP.proxy.final"
    else
      # If not biallelic, concatenate ":multi" to the 6th field and echo the modified line to the final file
      modified_line="${line//$field6/$field6:multi}"
      echo "$modified_line" >> "${SNAPTMP}/SNAP.proxy.final"
    fi
  done < "$file"
done


	awk 'FNR==NR { m[$3] = m[$3] " " $6 " " $7 " " $8; next } $2 == "." && ($1 in m) { print $1 m[$1], "proxy"; c += 1; next } { print } END { print "proxy: Matched", c, "SNP identifiers by proxy." > "/dev/stderr" }' "${SNAPTMP}/SNAP.proxy.final" ${SNAPTMP}/SNAP.TMP.proxy


}

update () {
#   Purpose: It is possible a query SNP went unmatched because it's RS_id is not up to date. This function will update the unmatched query SNPs to the latest RS_id for that SNP.
#            This output is piped back into the match() function to see if new matches can be made with the correct RS_id's.
#   Output: 
#          
#            Assuming the output piped to update() was:
#            
#            rs74196558 . . . .
#            rs59226740 . . . .
#            rs39372 . . . .
#            rs928473 . . . .
#           
#            update would output:
#
#            rs12895622 . . . .
#            rs12117927 . . . .
#            rs39372 . . . .
#            rs928473 . . . .


    cat </dev/stdin > ${SNAPTMP}/SNAP.TMP.update #Store output from position() or match() into the file SNAP.TMP.update
    echo -n "update: Updating $(awk '$2=="." { print $1 }' ${SNAPTMP}/SNAP.TMP.update | wc -l) SNPs to the latest dbSNP identifiers, please wait " > /dev/stderr
    for rs in $(awk '$2=="." { print $1 }' ${SNAPTMP}/SNAP.TMP.update) #access RS_id's of SNPs that are unmatched and iteratively update them.
    do
	rsid=$(echo $rs | perl -p -e "s/rs//g;") #sets the variable 'rsid' to be the RS_id withut the "rs"

	#curl -X GET --header "Accept: application/json" "https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/${rsid}" > ${SNAPTMP}/${rs}.dbsnp.json 2> ${SNAPTMP}/${rs}.dbsnp.log
	curl -X GET --header "Accept: application/json" "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/${rsid}" > ${SNAPTMP}/${rs}.dbsnp.json 2> ${SNAPTMP}/${rs}.dbsnp.log #removed beta from link
	#This retrieves the JSON file for the SNP and stores it in the ${rs}.dbsnp.log file.
	sleep 0.5 #delay before querying the server again.
	echo -n "." > /dev/stderr
    done
    echo "" > /dev/stderr
    grep -h merged_into ${SNAPTMP}/*.dbsnp.json | perl -p -e "s/\{\"refsnp_id\"\:\"([0-9]+).*merged_into\"\:\[\"([0-9]+).*/rs\1 rs\2/g;" > ${OUTPUT_FILE}.updated.txt #makes a file of the format "OLD_RSID NEW_RSID" 
    #if there is new version
    awk 'FNR==NR { m[$1]=$2; next } $1 in m { $1=m[$1]; c+=1 } { print } END { print "update: Updated",c,"SNP identifiers to most recent dbSNP identifiers, conversions written to updated.txt." > "/dev/stderr" }' ${OUTPUT_FILE}.updated.txt ${SNAPTMP}/SNAP.TMP.update #Produces the output with the new RS_ids to be piped back into match() to see if new matches can be made now.
}

clean () {
    rm ${SNAPTMP}/*
}




UCSC_DBSNP=/project/richards/ethan.kreuzer/snp156.db #this needs to chnage this back to $(dirname -- $0)/snp150.db format
SNAPTMP=SNAPTMP
mkdir -p ${SNAPTMP} #make new directory

#echo "Welcome to snappy! Select ancestries by typing its corresponding code in the second column.\n" > /dev/stderr


QSNP_F=$1
GWAS_F=$2
PLINK_REF_PANEL=$3
ANCESTRY_SELECTION=$4
PLINK_WINDOW=$5
PLINK_MIN_R2=$6
NUM_PROXYS=$7
OUTPUT_FILE=$8

pop_selection

init ${QSNP_F} ${GWAS_F} | match | update | match  | position | proxy ${PLINK_REF_PANEL} ${PLINK_WINDOW} ${PLINK_MIN_R2} ${NUM_PROXYS} > ${OUTPUT_FILE}

rm -r ${SNAPTMP}
