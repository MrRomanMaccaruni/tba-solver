// Stefano Sinigardi, Alessandro Fabbri // BSD License

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

bool Belongs_to(char c, std::string s){
  for (int i = 0; i < s.size(); i++){ if (c == s.at(i)) return true; }
  return false;
}

std::vector< std::vector<std::string> > Parse_file(std::ifstream &file_to_parse){
  // Parsing parameter to EDIT
	std::vector<std::string> fields({"number_of_steps","x_max","r","kernel_filename"});
	std::string separators="\t =,";
	std::string comment="#!;";
  // Internal variables
	std::string line;
	std::vector<std::string> tokens;
    std::vector< std::vector<std::string> > parsed;
  // Parsing loop
	while (!file_to_parse.eof()){
	  line.clear(); tokens.clear();
	  getline(file_to_parse, line);
	  boost::algorithm::trim(line);                                                                                   // remove leading/trailing spaces
	  if( !line.size() ) continue;                                                                                    // ignore empty lines
	  boost::algorithm::split(tokens, line, boost::algorithm::is_any_of(separators), boost::token_compress_on);       // split line to tokens, BOOST required
	  std::transform(tokens[0].begin(), tokens[0].end(), tokens[0].begin(), ::tolower);                               // convert to lower, NOT multi-byte char compatible
	  for(int i=0; i<tokens.size(); i++){                                                                             // remove inline comments
	  	if( Belongs_to(tokens[i][0],comment) ){ tokens.erase( tokens.begin()+i, tokens.end()); }
	  }
	  if( tokens.size() ){                                                                                            // ignore empty tokens
	  	for(int i=0; i<fields.size(); i++){
	  		if( tokens[0] == fields[i] ) {
	  			parsed.push_back( tokens );
	  		}
	  	}
	  }
	}
	return parsed;
}
