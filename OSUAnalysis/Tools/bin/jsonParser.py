import os
import json
import re
import subprocess

jsonName = 'Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON'
#jsonName = 'Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v2'

input  = open(jsonName+'.txt','r')
output = open(jsonName+'.C','w')

code = "bool isGood(int run, int lumi){\n"

for line in input.readlines():
   cert_list = json.loads(line)
   for run in sorted(cert_list):
      code += "if(run==" + run + " && ("
      for lumi_range in cert_list[run]:
         if( code[len(code)-1]==')' ):
             code += " ||";
         code += " (lumi>=" + str(lumi_range[0]) + " && lumi<=" + str(lumi_range[1]) + ")"
      code += ") ) return true;\n"

code += "return false; }"

luminosity_log = subprocess.Popen(['lumiCalc.py', '-c', 'frontier://LumiCalc/CMS_LUMI_PROD', '-i', jsonName+'.txt', 'overview', '-b', 'stable'], shell=False, stdout=subprocess.PIPE)

#luminosity_log = subprocess.Popen("""
# echo  '| 166860 |           36 |       804824.563 |               [1-21] |      389767.675 |' &&
# echo  '| 166861 |           50 |      1100317.288 |        [1-6], [8-13] |      176324.613 |' &&
# echo  ' ==  =  Total : ' &&
# echo  '| Delivered LS | Delivered (/mkb) | Selected LS | Recorded (/mkb) |' &&
# echo  '-------------------------------------------------------------------' &&
# echo  '|        57635 |    772727738.518 |       52540 |   702347412.382 |' &&
# echo  &&
# echo  ['WARNING] : run 166346 : total non-GOOD LS: UNKNOWN 167, BAD 0, SUSPECT 0' &&
# echo  ['WARNING] : run 165120 : total non-GOOD LS: UNKNOWN 16, BAD 0, SUSPECT 0'  &&
# echo  ['WARNING] : run 163332 : total non-GOOD LS: UNKNOWN 651, BAD 0, SUSPECT 0'
#""", shell=True, stdout=subprocess.PIPE)

code += "\n" 

counter = 0

for line in re.split('\n',luminosity_log.communicate()[0]):
   if re.search(' ==  =  Total :',line) != None :
      counter = 1
   if counter == 4 :
      code += "double luminosity(void){ return " + re.split('\s*\|\s*',line)[4] + ";}\n"
      break
   if counter!=0 :
      counter+=1

output.write(code);
output.close();

input.close()
