#!/usr/bin/python2 
'''
 This script does the work of creating single precision header and source files
 from their double precision version.

 There are three categories of files:
 1. All files (see list "fileNames")
 2. files where we can change all 'd_complex' with 'f_complex', 
 and also function names without fears ("filesALLDtoF"). Not all files are such.
 3. Files in single precision which have been modified by hand and must 
 not be overwritten  ("filesNoOverwrite")
 
 The conversion is performed in this way:
 1. Types are converted through a dictionary (for files where this is allowed)
 2. function names are found by searching all files for function declarations,
    and then replaced by the same name with "_f" appended
 3. global variable names are found by searching all files for their declaration,
    and then replaced by the same name with "_f" appended
 4. #define'd Numerical constant must also be converted to float,
    (for files where this is allowed)
 5. Library functions must also be converted to their float version 
    (for files where this is allowed)
 6. New filenames are the old filenames with 'sp_' prepended.
 7. sp_struct_c_def.h is created from struct_c_def.h, defining "f_complex"
 8. header guards are also changed, since a file with a new name is created.
    NOTICE that header guards must satisfy some criteria in order to be matched:

       headerGuard = os.path.basename(fileName).upper().replace('.','_')

'''
from sys import argv,exit
import os 
import re

# in autoMode, the user is not asked for permissions
# the script does not overwrite not-automatically written files
autoMode = False
if 'autoMode' in argv :
    autoMode = True
    argv.remove('autoMode')


# the script looks at output file timestamps, and by default does not 
# rewrite files
checkEverythingAnyway = False

if 'recheck' in argv :
    checkEverythingAnyway = True
    argv.remove('recheck')

silentMode = False
if 'silentMode' in argv :
    silentMode = True
    argv.remove('silentMode')

fileNames = [\
'OpenAcc/alloc_vars.c',\
'OpenAcc/alloc_vars.h',\
'OpenAcc/backfield.c',\
'OpenAcc/backfield.h',\
'OpenAcc/cayley_hamilton.h',\
'OpenAcc/fermion_force.c',\
'OpenAcc/fermion_force.h',\
'OpenAcc/topological_action.c',\
'OpenAcc/topological_action.h',\
'OpenAcc/topological_force.c',\
'OpenAcc/topological_force.h',\
'OpenAcc/fermion_force_utilities.c',\
'OpenAcc/fermion_force_utilities.h',\
'OpenAcc/fermionic_utilities.c',\
'OpenAcc/fermionic_utilities.h',\
'OpenAcc/fermion_matrix.c',\
'OpenAcc/fermion_matrix.h',\
'OpenAcc/inverter_full.c',\
'OpenAcc/inverter_full.h',\
'OpenAcc/inverter_multishift_full.c',\
'OpenAcc/inverter_multishift_full.h',\
'OpenAcc/ipdot_gauge.c',\
'OpenAcc/ipdot_gauge.h',\
'OpenAcc/matvecmul.h',\
'OpenAcc/md_integrator.c',\
'OpenAcc/md_integrator.h',\
'OpenAcc/plaquettes.c',\
'OpenAcc/plaquettes.h',\
'OpenAcc/rettangoli.c',\
'OpenAcc/rettangoli.h',\
'OpenAcc/single_types.h',\
'OpenAcc/stouting.c',\
'OpenAcc/stouting.h',\
'OpenAcc/struct_c_def.h',\
'OpenAcc/su3_measurements.c',\
'OpenAcc/su3_measurements.h',\
'OpenAcc/su3_utilities.c',\
'OpenAcc/su3_utilities.h',\
'DbgTools/dbgtools.h',\
'DbgTools/dbgtools.c',\
'Mpi/communications.h',\
'Mpi/communications.c',\
'Meas/gauge_meas.h',\
'Meas/gauge_meas.c',\
]

# files where we can change all 'd_complex' with 'f_complex'
# and also function names
# without fears
filesALLDtoF=[\
'OpenAcc/alloc_vars.c',\
'OpenAcc/alloc_vars.h',\
'OpenAcc/backfield.c',\
'OpenAcc/backfield.h',\
'OpenAcc/cayley_hamilton.h',\
'OpenAcc/fermion_force.c',\
'OpenAcc/fermion_force.h',\
'OpenAcc/fermion_force_utilities.c',\
'OpenAcc/fermion_force_utilities.h',\
'OpenAcc/topological_action.c',\
'OpenAcc/topological_action.h',\
'OpenAcc/topological_force.c',\
'OpenAcc/topological_force.h',\
'OpenAcc/fermion_matrix.c',\
'OpenAcc/fermion_matrix.h',\
'OpenAcc/ipdot_gauge.c',\
'OpenAcc/ipdot_gauge.h',\
'OpenAcc/matvecmul.h',\
'OpenAcc/plaquettes.c',\
'OpenAcc/plaquettes.h',\
'OpenAcc/rettangoli.c',\
'OpenAcc/rettangoli.h',\
'OpenAcc/single_types.h',\
'OpenAcc/stouting.c',\
'OpenAcc/stouting.h',\
'OpenAcc/struct_c_def.h',\
'OpenAcc/su3_measurements.c',\
'OpenAcc/su3_measurements.h',\
'OpenAcc/su3_utilities.c',\
'OpenAcc/su3_utilities.h',
'Mpi/communications.h',\
'Mpi/communications.c',\
'DbgTools/dbgtools.h',\
'DbgTools/dbgtools.c',\
'Meas/gauge_meas.h',\
'Meas/gauge_meas.c',\
]

# protected file lists
# these files will not be automatically overwritten
# confirmation will alway be asked singularly
# Possible reason(s) for including a file in this list:
# 1. the double precision file has modifications that 
#    must not be included in the single precision one
filesNoOverwrite=[\
'OpenAcc/fermion_force.c',\
]


dpFunctionNames = [] # new function names will just be dp function names + '_f' at the end
dpVariableNames = [] # new function names will just be dp function names + '_f' at the end
# as in struct_c_def.h
dpTypes = ['double_soa','dcomplex_soa','vec3_soa','vec3',\
        'su3_soa','thmat_soa','tamat_soa',\
        'global_vec3_soa','global_su3_soa',\
        'global_tamat_soa','global_thmat_soa',\
        'global_dcomplex_soa','global_double_soa',\
        'single_su3' , 'single_tamat', 'single_thmat']
# corresponding types in in sp_struct_c_def.h
spTypes = ['float_soa','fcomplex_soa','vec3_soa_f','vec3_f',\
        'su3_soa_f','thmat_soa_f','tamat_soa_f',\
        'global_vec3_soa_f','global_su3_soa_f',\
        'global_tamat_soa_f','global_thmat_soa_f',\
        'global_dcomplex_soa_f','global_double_soa_f',\
        'single_su3_f' , 'single_tamat_f', 'single_thmat_f']

# as in struct_c_def.h
dpTypes_t = [ dpType + '_t' for dpType in dpTypes]
spTypes_t = [ spType + '_t' for spType in spTypes]
dpTypes += dpTypes_t
spTypes += spTypes_t


dpToSpDict = dict(zip(dpTypes,spTypes))

# this function is needed to populate the list 'dpFunctIOnNames'


def findFunctionNames(lineRaw):
    returnTypes=['int','double','void','d_complex','vec3']
    foundFunction = False
    line = lineRaw.strip()
    for returnType in returnTypes:
        # find the name of the function
        reToMatch = '(?<='+returnType+'[ \*])' # match if preceded by returnType and ' ' or *
                                          # this re construct requires fixed length
        reToMatch += '[\s\*]*' # any space or asterisk or nothing (this is in the string )
        reToMatch += '\w+'   # any number >= 1 of characters in azAZ09_
        reToMatch += '(?=[\s]*[\(])' # match if followed by '('
        foundSomething = re.search(reToMatch,line)
        if foundSomething :
            rawFunctionName = line[foundSomething.start():foundSomething.end()]
            # the function name should be at the end of this
            # this second step is necessary because of the limitations of '(?<=...)'
            # since we can't match 'prefixes' of unknown length
            funcNameLocation = re.search('\w+$',rawFunctionName)
            # $ = end of rawFunctionName    
            dpFunctionName = rawFunctionName[funcNameLocation.start():funcNameLocation.end()]
            if dpFunctionName not in dpFunctionNames and dpFunctionName != 'main':
                if not silentMode :
                    print "Found function " + dpFunctionName
                dpFunctionNames.append(dpFunctionName)
            break



def findFirstVarName(text): # note: does not work with,e.g. 'int a,b;':
                            # only ONE var per declaration
    upperBoundaries = []
    for soaType in dpTypes:
        reToMatch = '(?<='+soaType+'[ \*])' # match if preceded by returnType and ' ' or *
                                          # this re construct requires fixed length
        reToMatch += '[\s\*]*' # any space or asterisk or nothing (this is in the string )
        reToMatch += '\w[\s\w]*'   # any number >= 1 of characters in azAZ09_
        reToMatch += '(?=[\s]*[;])' # match if followed by ';'
        foundSomething = re.search(reToMatch,text)
        if foundSomething:
            rawVarName = text[foundSomething.start():foundSomething.end()]
            upperBoundaries.append( foundSomething.end())
            # the function name should be at the end of this
            # this second step is necessary because of the limitations of '(?<=...)'
            # since we can't match 'prefixes' of unknown length
            # this should also remove 'const', '__restrict' i tak dalej
            varNameLocation = re.search('\w\w*$',rawVarName)
            # $ = end of rawFunctionName    
            dpVarName = rawVarName[varNameLocation.start():varNameLocation.end()]
            if dpVarName not in dpVariableNames: 
                if not silentMode:
                    print "Found Variable (", soaType, ')', dpVarName
                dpVariableNames.append(dpVarName)

    if len(upperBoundaries) == 0:
        return -1
    else:
        return min(upperBoundaries)



# collecting all function names
for fileName in fileNames:
    if not silentMode:
        print "----------\nChecking file" , fileName , '\n-----------'
    f  =open(fileName)
    fileLines = f.readlines()
    for line in fileLines:
        findFunctionNames(line)
    f.close()


dpFunctionNames.sort(key = len, reverse = True ) # here reverse is CRUCIAL 
                                                 # because some function names 
                                                 # contain other function names

functionNamesFile = open('functions_found.txt','w')
for foundFunction in dpFunctionNames:
    functionNamesFile.write(foundFunction + '\n')
functionNamesFile.close()
if not silentMode:
    print "Total functions found: ", len(dpFunctionNames)

# collecting all global variable names
if not silentMode:
    print "\n\n-----------------\nLooking for variable names\n-----------------\n"
f = open('OpenAcc/alloc_vars.c')
text = f.read()
foundSomething = True
index = 0
while foundSomething:
     res = findFirstVarName(text[index:])
     if res == -1:
         foundSomething = False
     else:
         index += res
f.close()


dpVariableNames.sort(key = len , reverse = True) # CRUCIAL
globalVarNamesFile = open('glvar_found.txt','w')
for foundGlobalVar in dpVariableNames:
    globalVarNamesFile.write(foundGlobalVar + '\n')
globalVarNamesFile.close()

dpVariableNames.append('phases') # for U1 used in dirac matrix
dpVariableNames.sort(key = len , reverse = True) # CRUCIAL

if not silentMode:
    print "Total global variables found: ", len(dpVariableNames)




# collected all double precision function names, now just replacing:
# - all function names
# - all relevant types

ans = '' # answer for 'overwrite file?'

if autoMode:
    ans = 'a'


fileNamesToChange = list(set(fileNames)-set(filesNoOverwrite))
fileNamesToChange += filesNoOverwrite

changedFiles = []

fileArgs = [ fileName.replace('sp_','') for fileName in argv ] 

for fileName in fileNamesToChange:
    # notice: if autoMode was in argv, it was previously removed 
    # so the case len(argv)>1 is the case when there are filenames passed as arguments 
    # and the case len(argv)==1 is the one without filenames passed as arguments
    if (len(argv) > 1 and fileName in fileArgs) \
            or (len(argv)==1 and fileName not in filesNoOverwrite):
        
        newFileName = os.path.dirname(fileName)+'/sp_'+os.path.basename(fileName)
        writeIt = True

        if fileName in filesNoOverwrite and autoMode:
            writeIt = False
        elif os.path.exists(newFileName):
            doubleFileModTime = os.path.getmtime(fileName)
            singleFileModTime = os.path.getmtime(newFileName)
            if singleFileModTime > doubleFileModTime and not checkEverythingAnyway:
                if silentMode:
                    print newFileName, "is already ok." 
                else:
                    print "File ", newFileName, " is newer than ",fileName," and won't be touched."
                
                writeIt = False

        elif os.path.exists(newFileName) and ( ans != 'a' or fileName in filesNoOverwrite):
            ans = ''
            while ans not in ['y','n','a']:
                if fileName in filesNoOverwrite:
                    print "\nWARNING: File " + fileName + " is in the 'protected' file list!"

                print "Overwrite file \'"+newFileName+"\'? (y=yes,n=no,a=yes to all)"
                ans = raw_input().lower()
                if ans == 'n':
                    writeIt = False



        if writeIt:
        
            f = open(fileName,'r')
            text = f.read()
            newText = str(text)
    
            # replacing function names
            for dpFunctionName in dpFunctionNames:
                reToMatch = '((?<=\W)|^)' # either preceded by the beginning of the string or a 
                                          # non-alphanumeric character
                reToMatch += dpFunctionName 
                reToMatch += '(?=(\W|$))' # either followed by the end of the string or a  
                                          # non-alphanumeric character
                newText = re.subn(reToMatch, dpFunctionName + '_f', newText)[0]
            # changing relevant (soa-like, arrays) types
            # note : this step could also change function names
            for dpType in dpToSpDict:
                reToMatch = '((?<=\W)|^)' # either preceded by the beginning of the string or a 
                                          # non-alphanumeric character
                reToMatch += dpType 
                reToMatch += '(?=(\W|$))' # either followed by the end of the string or a  
                                          # non-alphanumeric character
    
                newText = re.subn(reToMatch,dpToSpDict[dpType], newText)[0]
    
            # changing filenames in '#includes'
            for fileName2 in fileNames:
                reToMatch = '((?<=\W)|^)' # either preceded by the beginning of the string or a 
                                          # non-alphanumeric character
                reToMatch += os.path.basename(fileName2) 
                reToMatch += '(?=(\W|$))' # either followed by the end of the string or a  
                                          # non-alphanumeric character
                newText = re.subn(reToMatch,'sp_'+os.path.basename(fileName2),newText)[0]
            #taking care of header guards
            headerGuard = os.path.basename(fileName).upper().replace('.','_')
            if headerGuard not in newText:
                print "Warning, header guard \'" + headerGuard + "\' not found in file " + fileName
            newText = newText.replace(headerGuard, "SP_"+headerGuard)
            # taking care of global variables
            for dpVariableName in dpVariableNames:
                reToMatch = '((?<=\W)|^)' # either preceded by the beginning of the string or a 
                                          # non-alphanumeric character
                reToMatch += dpVariableName
                reToMatch += '(?=(\W|$))' # either followed by the end of the string or a  
                                          # non-alphanumeric character
    
                newText = re.subn(reToMatch,dpVariableName + '_f', newText)[0]
    
    
    
        
        
        
            # it may happen that two transformations appear on the same symbol, 
            # and an '_f_f' is appended instead of just '_f'
            #newText = newText.replace('_f_f','_f')
            #newText = newText.replace('_f_f','_f')
        
        
            newText = newText.replace('deltas_Omelyan','deltas_Omelyan_f')
            newText = newText.replace('DOUBLE PRECISION VERSION','SINGLE PRECISION VERSION')
            allSubst = []
            allSubst.append(('d_complex','f_complex'))
            allSubst.append(('MPI_DOUBLE','MPI_FLOAT'))
            allSubst.append(('double','float'))
            allSubst.append(('conj','conjf'))
            allSubst.append(('creal','crealf'))
            allSubst.append(('cimag','cimagf'))
            allSubst.append(('cos','cosf'))
            allSubst.append(('sin','sinf'))
            allSubst.append(('exp','expf'))
            allSubst.append(('sqrt','sqrtf'))
            allSubst.append(('pow','powf'))
            allSubst.append(('fabs','fabsf'))
            allSubst.append(('C_ZERO','C_ZEROF'))
            allSubst.append(('C_ONE','C_ONEF'))
            allSubst.append(('RHO','RHOF'))
            allSubst.append(('ONE_BY_THREE','ONE_BY_THREEF'))
            allSubst.append(('ONE_BY_SIX','ONE_BY_SIXF'))
            allSubst = dict(allSubst)
    
    
            if fileName in filesALLDtoF:
                newText = newText.replace('%lf','%f')
                newText = newText.replace('%.18lf','%f')
                for subst in allSubst:
                    reToMatch = '((?<=\W)|^)' # either preceded by the beginning of the string or a 
                                              # non-alphanumeric character
                    reToMatch += subst
                    reToMatch += '(?=(\W|$))' # either followed by the end of the string or a  
                                          # non-alphanumeric character
    
                    newText = re.subn(reToMatch,allSubst[subst], newText)[0]
    
                # find number to convert to float
                reToMatch = '[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9]+)?'
                reNotToMatch = '%[-+]?[0-9]*\.[0-9]+([eE][-+]?[0-9]+)?' # no match for format
                                                                        # specifiers 
                listend = [0]  # end  f each match
                for m in re.finditer(reToMatch,newText):
                    listend.append(m.end())
                for m in re.finditer(reNotToMatch,newText):
                    if m.end() in listend:
                        listend.remove(m.end())
                newText2 = ''
           
                for i in range(len(listend)-1):
                    newText2 +=newText[listend[i]:listend[i+1]] + 'f'
                newText2 += newText[listend[-1]:]
                newText = newText2
    
    
    
           
            
            # adding 'typedef double complex d_complex' in sp_struct_c_def.h
            if  'struct_c_def.h' in fileName:
                if not silentMode:
                    print "Adding \'typedef float complex f_complex\' to ", newFileName
                newText = newText.replace('//TYPEDEF_FLOAT_COMPLEX','typedef float complex f_complex;' )
        
            if os.path.exists(newFileName):

                if not silentMode:
                    print "Checking if ", newFileName, " must be modified..."
                oldFile = open(newFileName,'r')
                oldText = oldFile.read()
                oldFile.close()
                if newText != oldText:
                    if not silentMode:
                        print "... file must be modified. Writing file ", newFileName , " ..."
                else:
                    if not silentMode:
                        print "File ", newFileName, " has no changes and won't be touched."
                    writeIt = False
                
                
                









            if writeIt:
                if not silentMode:
                    print "Writing file ", newFileName, "..."
                newFile = open(newFileName,'w')
                newFile.write(newText)
                newFile.close()
                changedFiles.append(newFileName)
        f.close()



if not silentMode :     
    print "\n\nRESUME: Changed files:"
    for changedFile in changedFiles:
        print changedFile
    if len(changedFiles) == 0 :
        print "None!"
    
else :
    strToPrint = "Changed files:"
    for changedFile in changedFiles:
        strToPrint += "\n" + changedFile
    if len(changedFiles) == 0 :
        strToPrint += " None!"
    print strToPrint
 
