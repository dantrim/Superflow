#pragma once

#include <string>
#include <sstream>
#include <vector>

namespace sflow {

    std::string pad_width(std::string s, int width)
    {
        std::stringstream stream_;
        int pad = width - s.size();
        stream_ << s;
        while (pad > 0) {
            stream_ << " ";
            pad--;
        }
        return stream_.str();
    }


    void split(const std::string &s, char delim, std::vector<std::string> &elems) // harder to use
    {
        std::stringstream ss(s);
        std::string item;
        while (getline(ss, item, delim)) {
            elems.push_back(item);
        }
    }

    std::vector<std::string> split(const std::string &s, char delim) // use this form
    {
        std::vector<std::string> elems;
        split(s, delim, elems);
        return elems;
    }

    std::string trim(std::string s)
    {
        std::string val = s.erase(s.find_last_not_of(" \n\r\t") + 1);
        val.erase(0, s.find_first_not_of(" \n\r\t"));
        return val;
    }

};
