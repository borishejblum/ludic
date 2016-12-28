#include <Rcpp.h>
using namespace Rcpp;

//' Splitting a character string in C++
//' 
//'@rdname strsplitC
//'
//'@description \code{strsplitC_memoryPbs} can have memory issues...
//' 
//'@examples
//'strsplitC_memoryPbs(c(";aaa;bb ; cccc;ee;"), sep=";")
//'@export
// [[Rcpp::export]]
std::vector<std::string> strsplitC_memoryPbs(std::string s, char sep){
  
  int len = s.size();
  int nb_split = 0;
  IntegerVector split_index;
  split_index[0] = -1;
  std::vector<std::string> out;
  for(int i=0; i < len; i++ ){
    if(s[i] == sep){
      nb_split += 1;
      split_index[nb_split]=i;
      if(i!=0){
        out.push_back(s.substr(split_index[nb_split-1]+1, i-(split_index[nb_split-1]+1)));
      }
    }
  }
  if(split_index[nb_split] != len-1){
    out.push_back(s.substr(split_index[nb_split]+1, len));
  }
  return(out);
}

//'@rdname strsplitC 
//'@description \code{strsplitC_safe} implementation is safer in its memory treatment
//' 
//'@examples
//'strsplitC_safe(c(";aaa;bb ; cccc;ee;"), sep=";")
//'strsplitC_safe(c(""), sep=";")
//'strsplitC_safe(c("a"), sep=";")
//'strsplitC_safe(c("bb ; ee"), sep=";")
//'strsplitC_safe(c("bb ee"), sep=";")
//'@export
// [[Rcpp::export]]
std::vector<std::string> strsplitC_safe(std::string s, char sep){
  
  int len = s.size();
  int nb_split = 0;
  
  for(int i=0; i < len; i++ ){
    if(s[i] == sep){
      nb_split += 1;
    }
  }
  
  
  std::vector<std::string> s_splitted = std::vector<std::string>(nb_split+1);
  int which_split = 0;
  IntegerVector split_index = IntegerVector(nb_split+1);
  split_index[0] = -1;
  int nb_empty =0;
  
  for(int i=0; i < len; i++ ){
    if(s[i] == sep){
      which_split += 1;
      split_index[which_split] = i;
      s_splitted[which_split-1] = (s.substr(split_index[which_split-1]+1, i-(split_index[which_split-1]+1)));
      if(s_splitted[which_split-1] == ""){
        nb_empty += 1;
      }
    }
  }
  
  
  if(split_index[which_split] != len-1){
    s_splitted[which_split] = s.substr(split_index[which_split]+1, len);
  }
  if(s_splitted[which_split] == ""){
    nb_empty += 1;
  }
  
  if(nb_empty>0){
    int which_notempty = 0;
    std::vector<std::string> out = std::vector<std::string>(nb_split + 1 - nb_empty);
    for(int i=0; i<nb_split + 1 ; i++){
      if(s_splitted[i] != ""){
        which_notempty += 1;
        out[which_notempty-1] = s_splitted[i];
      }
    }
    return(out);
  }
  else{
    return(s_splitted);
  }
}


//' Splitting a character string in C++
//' 
//'Splitting a subset of each character string from character 
//'number \code{i} to character \code{l} on characters \code{";"}.
//' 
//'@param Adates a matrix of character strings to be subsetted and splitted
//'
//'@param i an integer indicating the begining of the subset of the strings
//'
//'@param l an integer indicating the end of the subset of the strings
//' 
//'@examples
//'stringdate_test(matrix(c("20NOV2008:21NOV2008", 
//'                         "20NOV2008:21NOV2008 ; 10APR2008:12APR2008", "", ""), 
//'                       ncol=2, nrow=2), 1,0)
//'
//'@keywords internal                          
//'                          
//'@export
// [[Rcpp::export]]
int stringdate_test(StringMatrix Adates, int i, int l){
  std::string datesA_temp_string = Rcpp::as<std::string>(Adates(i, l));
  std::vector<std::string> datesA_stringvect = strsplitC_memoryPbs(datesA_temp_string, ';' );
  int nb_datesA = datesA_stringvect.size();
  return(nb_datesA);
}

//' C++ function for comparing dates encoded as character stings
//' 
//'@param s1 a date as a character string
//'
//'@param s2 a second date as a character string
//'
//'@param fmt1 a fomat to read the date \code{s1}
//'
//'@param fmt2 a fomat to read the date \code{s2}
//' 
//'@return \code{TRUE} if \code{s1} precedes \code{s2}, \code{FALSE} otherwise.
//' 
//' @examples
//'dateread("20NOV2008 ", "%d%B%Y")
//'dateread(" JAN 10 2008", "%B %d% Y")
//'dateread(" Jan 10 2008 12:00AM", " %B  %d %Y")>=dateread(" Jan 26 2008 12:00AM", " %B  %d %Y")
//'datetest(" Nov 3 2008 12:00AM", "04NOV2008", " %B %d %Y", "%d%B%Y")
//'datetest(" Nov 5 2008 12:00AM", "04NOV2008", " %B %d %Y", "%d%B%Y")
//'
//'
//'@keywords internal
//'
//'@export
// [[Rcpp::export]]
bool datetest(String s1, String s2, std::string fmt1, std::string fmt2){
  //Rcout << (Date(s1, fmt1) < Date(s2, fmt2));
  bool out = (Date(s1, fmt1) < Date(s2, fmt2)); 
  return(out);
}

//' C++ function for reading dates from character stings
//' 
//'@param s a date as a character string
//'
//'@param fmt a fomat to read the date
//'
//'@examples
//'dateread(" Nov 10 2008 12:00AM", "%B%d%Y")
//'dateread(" Nov 10 2008 12:00AM", " %B  %d %Y")
//'dateread(" Nov 10 2008 12:00AM", " %d %B %Y")
//'
//'@keywords internal
//'
//'@export
// [[Rcpp::export]]
Date dateread(String s, std::string fmt){
  return(wrap(Date(s, fmt)));
}


