#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <fstream>
#include <string>
#include <map>
#include <vector>
#include "tinyxml.h"

namespace levelset {
    
    class Value {
        
        std::string s;
        std::string desc;
        static int index_;
        int myindex;
        std::string units;
        
    public:
        
        Value(void) : s(""), myindex(index_++), units("") {}
        Value(const std::string q) : myindex(index_++), units("") {s = eatWhite(q);} 
        
        std::string operator=(const std::string q) {s = q; s = eatWhite(s); return s;}
        
        operator int() const {return atoi(s.c_str());} 
        operator double() const {return atof(s.c_str());} 
        operator char() const {return s[0];} 
        //   operator bool() const {return s == std::string("true") || s[0] == '1';} 
        operator std::string() const {return s;}
        
        bool isTrue(void) const {return s == std::string("true") || s[0] == '1';} 
        
        void SetDescription(const std::string q) {desc = q;} 
        void SetDescription(const char* c) {desc = std::string(c);} 
        std::string Description(void) {return desc;} 
        
        void SetUnits(const std::string u);
        std::string Units(void) {return units;}
        
        int Index(void) const {return myindex;} 
        
    private:
        
        std::string eatWhite(std::string q) const 
        {q.erase(0,q.find_first_not_of(" ")); return q;}
        
        friend class SortByIndex;
    };
    
    class IPElement {
        
        std::string name;
        std::map<std::string,Value> attrib;
        std::multimap<std::string,IPElement> child;
        
    public:
        
        IPElement(std::string nm) : name(nm) {}
        IPElement(void) {}
        void SetName(std::string nm) {name = nm;}
        std::string GetName(void) {return name;}
        void AddAttribute(const std::string label, const std::string value);
        Value GetAttribute(const std::string label);
        bool IsAttribute(const std::string label) const;
        Value GetValue(const std::string unit = "");
        std::multimap<std::string,IPElement> GetChildren(const std::string label = "!all");
        IPElement GetChild(const std::string label) const;
        unsigned int ChildCount(void) const {return (unsigned int)(child.size());}
        void AddChild(const IPElement& elem);
        friend std::ostream& operator<<(std::ostream& s, const IPElement& elem);
    };
    
    class InputParams {
        
        typedef std::map<std::string,Value> ipmap;
        ipmap values;
        ipmap updates;
        std::string fname;
        bool maketemplate;
        std::string bname;
        bool isXML;
        TiXmlDocument xmldoc;
        IPElement root;
        
    public:
        
        InputParams(char* filename);
        InputParams(int nargs, char* argv[]);
        
        template <class T>
        T      GetParam(const char* name, const T defalt, const char* desc = NULL);
        
        int    GetIntParam(const char* name);
        double GetDoubleParam(const char* name);
        double GetDoubleParamWithUnits(const char* name, const char* unit);
        void   GetCharParam(const char* name, char* ans);
        bool   GetBooleanParam(const char* name);
        int    GetIntParam(const char* name, const int defalt, const char* desc = NULL);
        double GetDoubleParam(const char* name, const double defalt,
                              const char* desc = NULL);
        void   GetCharParam(const char* name, char* ans, const char* defalt,
                            const char* desc = NULL);
        std::string GetStringParam(const char* name, const char* defalt = NULL, const char* desc = NULL);
        bool   GetBooleanParam(const char* name, const bool defalt, const char* desc=NULL);
        std::vector<IPElement> GetParamList(const char* name);
        IPElement GetElement(const char* name);
        IPElement GetElementTree(void) {return root;}
        
        void   CloseTemplate(void);
        void   AddParam(const char* buffer);
        
        char*  Basename(char* name) const {strcpy(name, bname.c_str()); return name;} 
        
        void WriteXML(const char* name = NULL);
        
    private:
        
        void ReadFile(const char* fname);   
        void ExtractNode(TiXmlNode* parent, IPElement& elem);
        Value GetParamValue(const char* name);
        Value GetParamValue(const char* name, const std::string defalt,const char* desc);
    };
    
}
#endif
