# CrazyBioUtils
These simple scripts can be combined into pipeline to address complex bioinformatics problems. The usage of each script can be displayed by typing ```xxx.pl -h```:
```shell
Start blast2go...
SYSNOPSIS
blast2go.pl [options] [options] -i example.fa -o example

 Options:
    -i --in_fa       sequences for blast2go in fa format
    -s --go-slim     converts the annotations into the GoSlim version
                     it may be wrong to in convert
    -r --iprscan     <path/to/adirectory/with/ipr_xml_files>
                     Note: not completed
    -p --progress    
    -o --out-folder  output the results to folder,default ./
    -b --b2g4pipe    <path/to/b2g4pipe/>, default:'~/blast/blast2go/b2g4pipe'
```
