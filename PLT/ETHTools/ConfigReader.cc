#include "ConfigReader.h"


ConfigReader::ConfigReader(const char* f, const char* d){
  ifstream fin;
  ofstream fout;

  strncpy(fileName,f,100);

  ifstream ftest(f, ios::in);// somehow ifstream doesn't work like I think it does
  if( !ftest ){
    cout << "unable to open config file " << f << endl;
	 if(strlen(d)>0){
		fin.open(d, ios::in);
		if( fin.is_open() ){
		  cout << "reading from " << d << endl;
		  fout.open(f);
		  fout << "#copied from default " <<  d << endl;
		}else{
		  cout << "running without config file" << endl;
		}
	 }
  }else{
	 fin.open(f, ios::in);
	 cout << "reading config "<< f << endl;
  }

  if (fin.is_open()) {
    string sLine;
    while( getline(fin,sLine) ){
		if (fout.is_open()) fout << sLine << endl;
      unsigned int i;
      ConfigItem a;

		a.line=sLine;
		a.id="void";
		a.idx1=kNoIndex;  // assume no index
		a.idx2=kNoIndex;


      vector<string> tokens;
      Tokenize(sLine, tokens,"[,]");

		if(tokens.size()>0){
		  //debugging
		  /*
		  cout << tokens.at(0);
		  for(int j=1; j<tokens.size(); j++){ cout << ":"<<tokens.at(j); }
		  cout << endl;
		  */
		  
		  a.id  =tokens.at(0);
		  if(tokens.size()>1){
			 i=1;
			 if((tokens.size()>3)&&(tokens.at(1)=="[")){
				a.idx1=atoi(tokens.at(2).c_str());
				i=3;
				if((tokens.size()>5)&&(tokens.at(3)=="'")){
				  a.idx2=atoi(tokens.at(4).c_str());
				  i=5;
				}
				if(tokens.at(i)!="]"){
				  cout << "error parsing line :" << sLine << endl 
						 << "closing ] expected instead of " << tokens.at(i) << endl;
				  break;
				}else{
				  i++;  // ok, proceed to the values.
				}
			 }
			 while(i<tokens.size()){
				a.values.push_back(tokens.at(i));
				i++;
			 }
		  }// size>1
		} // size >0
		// create one item per line, maybe void for comment lines
		items.push_back(a);
	 }// while getLine
  }// file open
  
  if(1){
    cout << "found " << items.size() << endl;
    for(unsigned int i=0; i<items.size(); i++){
      cout << i <<" " << items.at(i).id << " " << items.at(i).values.size() 
			  << " " << items.at(i).idx1 << ","<< items.at(i).idx2 << endl;
    }
  }
}


ConfigReader::~ConfigReader(){
  cout <<" bye bye "<< endl;
}

void ConfigReader::Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& symbols)
{
  // white spaces
  const string whitespace=" ";
  string delimiters=symbols+whitespace;

  // Skip whitespaces at beginning.
  string::size_type lastPos = str.find_first_not_of(whitespace, 0);
  // Find first delimiter
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
  {
  if(lastPos==pos){
  // found a single character symbol
  tokens.push_back(str.substr(lastPos, 1));
  pos++;  // move on to next character
      }else{
        // Found a string, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
      }
      // Skip whitespace.  Note the "not_of"
      lastPos = str.find_first_not_of(whitespace, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
  }
}


ConfigItem * ConfigReader::findItem(const char* name){
  for (unsigned int i=0; i<items.size(); i++){
    if(items.at(i).id==name){
      return &items.at(i);
    }
  }
  return 0;
}


ConfigItem * ConfigReader::findItem(const char* name, int idx){
  for (unsigned int i=0; i<items.size(); i++){
    if((items.at(i).id==name)&&(items.at(i).idx1==idx)){
      return &items.at(i);
    }
  }
  return 0;
}

void ConfigReader::warnLength(ConfigItem* a, int nRequest){
  cout << "ConfigReader> *** length mismatch: ";
  cout << "item "<< a->id ;
  if (!(a->idx1==kNoIndex)){
	 cout << "[" << a->idx1;
	 if (!(a->idx2==kNoIndex)){
		cout <<","<<a->idx2;
	 }
	 cout <<"]";
  }
  cout <<" has " <<a->values.size() << " values, " 
		 << nRequest << " were requested" << endl;
}


/***************************************************************************/

/* single integer: name value*/
void ConfigReader::get(const char* name, int &value, int defaultValue){
  ConfigItem * a= findItem(name);
  if( a ){
    value=atoi(a->values.at(0).c_str());
  }else{
    value=defaultValue;
  }
}

/* single indexed integer:  name[index] value*/
void ConfigReader::get(const char* name, const int idx, 
		      int &value, int defaultValue){
  ConfigItem * a= findItem(name,idx);
  if( a ){
    value=atoi(a->values.at(0).c_str());
  }else{
    value=defaultValue;
  }
}


/* single double:  name value*/
void ConfigReader::get(const char* name, double &value, double defaultValue){
  ConfigItem * a= findItem(name);
  if( a ){
    value=atof(a->values.at(0).c_str());
  }else{
    value=defaultValue;
  }
}

/* single indexed double:  name[index] value*/
void ConfigReader::get(const char* name, const int idx, 
		      double &value, double defaultValue){
  ConfigItem * a= findItem(name,idx);
  if( a ){
    value=atof(a->values.at(0).c_str());
  }else{
    value=defaultValue;
  }
}


/* array of integers:  name value_1 value_2 .... value_nval*/
void ConfigReader::geta(const char* name, const unsigned int nval, int* v){
  ConfigItem * a= findItem(name);
  if( a ){
    if(a->values.size()==nval){
      for(unsigned int i=0; i<nval; i++){
		  v[i]=atoi(a->values.at(i).c_str());
      }
    }else{ // length mismatch
		warnLength(a, 1);
    }
  }else{ // not found
    //    what?
  }
};

/* indexed array of integers:  name[index] value_1 value_2 .... value_nval*/
void ConfigReader::geta(const char* name, const int idx, const unsigned int nval, int* v){
  ConfigItem * a= findItem(name, idx);
  if( a ){
    if(a->values.size()==nval){
      for(unsigned int i=0; i<nval; i++){
	v[i]=atoi(a->values.at(i).c_str());
      }
    }else{ // length mismatch
		warnLength(a, nval);
    }
  }else{ // not found
    //    what?
  }
};

/* array of doubles:  name value_1 value_2 .... value_nval*/
void ConfigReader::geta(const char* name, const unsigned int nval, double* v){
  ConfigItem * a= findItem(name);
  if( a ){
    if(a->values.size()==nval){
      for(unsigned int i=0; i<nval; i++){
		  v[i]=atof(a->values.at(i).c_str());
      }
    }else{ // length mismatch
		warnLength(a, 1);
    }
  }else{ // not found
    //    what?
  }
};


/* indexed array of doubles:  name[index] value_1 value_2 .... value_nval*/
void ConfigReader::geta(const char* name, const int idx, const unsigned int nval, double* v){
  ConfigItem * a= findItem(name, idx);
  if( a ){
    if(a->values.size()==nval){
      for(unsigned int i=0; i<nval; i++){
		  v[i]=atof(a->values.at(i).c_str());
      }
    }else{ // length mismatch
		warnLength(a, nval);
    }
  }else{ // not found
    //    what?
  }
};






/* indexed array of doubles:  name[index] value_1 value_2 .... value_nval*/
void ConfigReader::updatea(const char* name, int idx, const unsigned int nval, 
			   double *v, const char* format){
  char buf[101];
  cout << "updating" << endl;
  ConfigItem * a=findItem(name,idx);
  if( a ){ 
    if(a->values.size()==nval){
      cout << "old:  " << a->line << endl;
      snprintf(buf,100,"%s[%d] ",name,idx);
      a->line=buf;
      for(unsigned int i=0; i<nval; i++){
	snprintf(buf,100,format,v[i]);
	a->values[i]=buf;
	a->line+=" ";
	a->line+=buf;
      }
      cout << "new:  " << a->line << endl;
    }
  }
}


void ConfigReader::updatea(const char* name, const unsigned int nval, 
			   double *v, const char* format){
  char buf[101];
  cout << "updating" << endl;
  ConfigItem * a=findItem(name);
  if( a ){ 
    if(a->values.size()==nval){
      cout << "old:  " << a->line << endl;
      snprintf(buf,100,"%s ",name);
      a->line=buf;
      for(unsigned int i=0; i<nval; i++){
	snprintf(buf,100,format,v[i]);
	a->values[i]=buf;
	a->line+=" ";
	a->line+=buf;
      }
      cout << "new:  " << a->line << endl;
    }
  }
}


/* indexed array of ints:  name[index] value_1 value_2 .... value_nval*/
void ConfigReader::updatea(const char* name, int idx, const unsigned int nval, 
			   int *v, const char* format){
  char buf[101];
  cout << "updating" << endl;
  ConfigItem * a=findItem(name,idx);
  if( a ) {
	 if(a->values.size()==nval){
		cout << "old:  " << a->line << endl;
		snprintf(buf,100,"%s[%d] ",name,idx);
		a->line=buf;
		for(unsigned int i=0; i<nval; i++){
		  snprintf(buf,100,format,v[i]);
		  a->values[i]=buf;
		  a->line+=" ";
		  a->line+=buf;
		}
		cout << "new:  " << a->line << endl;
	 }
  }else{
	 cout << name << "not found " << endl;
  }

}
/* array of ints:  name value_1 value_2 .... value_nval*/
void ConfigReader::updatea(const char* name, const unsigned int nval, 
			   int *v, const char* format){
  char buf[101];
  cout << "updating" << endl;
  ConfigItem * a=findItem(name);
  if( a ) {
    if(a->values.size()==nval){
      cout << "old:  " << a->line << endl;
      snprintf(buf,100,"%s ",name);
      a->line=buf;
      for(unsigned int i=0; i<nval; i++){
	snprintf(buf,100,format,v[i]);
	a->values[i]=buf;
	a->line+=" ";
	a->line+=buf;
      }
      cout << "new:  " << a->line << endl;
    }
  }else{
    cout << name << "not found " << endl;
  }
  
}





void ConfigReader::rewrite(){
  ofstream fout(fileName);
  
  for(unsigned int i=0; i<items.size(); i++){
	 cout << items[i].line << endl;
	 fout << items[i].line << endl;
  }
  fout.close();
  
}
