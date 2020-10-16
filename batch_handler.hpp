//
//  batch_handler.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 5/15/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#define batch_handler_hpp

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
using namespace std;
using namespace boost::filesystem;

vector<string> batch_handler(string dir){
    vector<string> files;
    vector<pair<string,int>> dirs;
    string csv = "csv",file;
    path p (dir), output_folder = p / "/output_folder_";
    if(is_regular_file(p)) 
        output_folder = p.parent_path() / "/output_folder_";
    try{
        if (exists(p)){
            int output_dir_num = 0;
            if(exists(output_folder) && is_directory(output_folder)){
                while(is_directory(output_folder)){
                    output_dir_num += 1;
                    string name = "/output_folder_" + to_string(output_dir_num);
                    if(!is_regular_file(p))output_folder = p / name;
                    else output_folder = p.parent_path() / name;
                }
            }
            create_directory(output_folder);
            files.push_back(output_folder.string());
            if (is_regular_file(p)){
                file = p.string();
                size_t csv_check = file.find(csv);
                if(csv_check != string::npos) 
                    files.push_back(file);
                return files;
            }else if (is_directory(p)){
                dirs.push_back(make_pair(p.string(),0));
                recursive_directory_iterator x(p), end;
                remove( x->path().string() + "/temp.csv" );
                while(x != end){
                    string filename = x->path().filename().string();
                    if (filename.find("output_") != string::npos ||filename.find(".") == 0||filename.find("temp_v_FOM") != string::npos)
                    {
                        x.no_push();
                        ++x;
                        continue;
                    }
                    if (is_regular_file(x->path())){
                        file = x->path().string();
                        size_t csv_check = file.find(csv);
                        if(csv_check != string::npos){
                            files.push_back(file);
                            dirs.back().second += 1;
                        }
                    }else if (is_directory(x->path())){
                        dirs.push_back(make_pair(x->path().string(),0));
                    }else
                        cout << p << " does not exist, try entering full path\n";
                    ++x;
                }
            }
            else cout << p << " exists, but is not a regular file or directory, try entering full path\n";
        }
        else
            cout << p << " does not exist, try entering full path\n";
    }
    catch (const filesystem_error& ex){
        cout << ex.what() << '\n';
    }
    cout<<"----- Batch Files Info -----"<<endl;
    cout<<"Directories found: "<<(dirs).size()<<endl;
    cout<<"CSV files found: "<<(files).size()-1<<endl;
    cout<<"Output directory created: "<<output_folder<<endl;
    cout<<"-- Directories Content --"<<endl;
    for(auto i = dirs.begin(); i != dirs.end(); ++i){
        cout<<i->first<<endl;
        cout<<"Contained "<< i->second <<" CSV files"<<endl<<endl;
    }
    cout<<"----------------------------"<<endl;
    return files;
}
