/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include <iostream>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <tuple>
#include <sstream>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

// #include "H5Cpp.h"
// using namespace H5;
#include "hdf5.h"

#include "libx/File.hpp"
#include "libx/GZFile.hpp"
#include "libx/String.hpp"
#include "libx/System.hpp"
#include "libx/Timer.hpp"

#include "CLI11.hpp"

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long uint64;

#define DEBUG

static constexpr int BIN_SIZE = 200;
static constexpr float SAMPLE_RATIOS[] = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
static constexpr float DATA_RATIO = 0.05;
static constexpr int MIN_BIN_THRE = 5;

// Store raw data from gene expression file
struct ExpLine
{
    ExpLine(uint32 cx, uint32 cy, uint32 gene, uint64 mid):
            cx(cx), cy(cy), gene(gene), mid(mid)
    {
    }
    uint32 cx;
    uint32 cy;
    uint32 gene;
    uint64 mid;
};

// Store raw data from tissue.gef
typedef struct exp_t
{
    uint32 x;
    uint32 y;
    uint32 cnt;
} exp_t;

typedef tuple<uint32, uint64> GeneMID;
typedef tuple<uint32, uint32, uint32, uint64> CoorGeneMID;

typedef unordered_map<uint64, map<GeneMID, uint32>> SatStatMap;
typedef unordered_map<uint64, set<CoorGeneMID>> MIDStatMap;

struct Arguments
{
    string expFile;
    string outDir;
    string tissueFile;

    vector<string> annoStatFiles;
    vector<string> bcMapStatFiles;

    string outSatFile;
};
Arguments arguments;

// Bin path, for executing plot.py
fs::path binPath;

static const std::string version = "2.3.5";

void save_error_code(int errcode, string info)
{
    std::ostringstream ostr;
    std::time_t        t =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    ostr << std::put_time(std::localtime(&t), "[%Y-%m-%d %H:%M:%S] SAW-A8000") << errcode
         << ": " << info << endl;

    // std::cerr<<ostr.str();

    string filename("errcode.log");
    libx::writeFile(ostr.str(), filename, true);
}

bool parseArgs(int argc, char** argv)
{
    libx::Timer timer;
    // Parse the command line parameters.
    CLI::App app{ "Sequencing Saturation and Plotting." };
    app.footer("saturation version: " + version);
    app.get_formatter()->column_width(40);

    // Required parameters
    app.add_option("-I,-i", arguments.expFile,
                   "Input gene expression matrix contains reads count information");
    app.add_option("-O,-o", arguments.outDir, "Output directory");
    app.add_option("--tissue", arguments.tissueFile, 
                   "Gene expression matrix only contians data under the tissue region");
    string annoStatFiles;
    app.add_option("--summary", annoStatFiles, 
        "Annotation summary file, support multi-files seperated by ','");
    string bcMapStatFiles;
    app.add_option("--bcstat", bcMapStatFiles, 
        "Barcode mapping statistics report, support multi-files seperated by ','");

    app.callback([&]() {
        if (arguments.expFile.empty())
        {
            save_error_code(1, "-I is missing.");
            exit(-1);
        }
        if (!fs::exists(arguments.expFile))
        {
            save_error_code(2, "not exists file: " + arguments.expFile);
            exit(-1);
        }
        if (arguments.outDir.empty())
        {
            save_error_code(1, "-O is missing.");
            exit(-1);
        }
        if (!fs::exists(arguments.outDir))
        {
            save_error_code(2, "not exists path: " + arguments.outDir);
            exit(-1);
        }
        if (arguments.tissueFile.empty())
        {
            save_error_code(1, "--tissue is missing.");
            exit(-1);
        }
        if (!fs::exists(arguments.tissueFile))
        {
            save_error_code(2, "not exists file: " + arguments.tissueFile);
            exit(-1);
        }
        if (annoStatFiles.empty())
        {
            save_error_code(1, "--summary is missing.");
            exit(-1);
        }
        if (bcMapStatFiles.empty())
        {
            save_error_code(1, "--bcstat is missing.");
            exit(-1);
        }
    });
        
    try 
    {
        app.parse(argc, argv);
    } 
    catch (const CLI::ParseError &e) 
    {
        save_error_code(1, "invalid parameters.");
        app.exit(e);
        return false;
    }

    arguments.outSatFile = arguments.outDir + "/sequence_saturation.tsv";
    // Check files is exist
    vector<string> temp;
    libx::split(annoStatFiles, ',', temp);
    for (auto& f : temp)
    {
        if (f.empty()) continue;
        if (!fs::exists(f))
        {
            save_error_code(2, "not exists file: " + f);
            exit(-1);
        }
        arguments.annoStatFiles.push_back(f);
    }
    temp.clear();
    libx::split(bcMapStatFiles, ',', temp);
    for (auto& f : temp)
    {
        if (f.empty()) continue;
        if (!fs::exists(f))
        {
            save_error_code(2, "not exists file: " + f);
            exit(-1);
        }
        arguments.bcMapStatFiles.push_back(f);
    }
    if (arguments.annoStatFiles.empty())
    {
        save_error_code(2, "--summary has no valid files.");
        exit(-1);
    }
    if (arguments.bcMapStatFiles.empty())
    {
        save_error_code(2, "--bcstat has no valid files.");
        exit(-1);
    }

    binPath = argv[0];
    binPath = fs::absolute(binPath).parent_path();

    return true;
}

void parseGefTissue(string tissueFile, unordered_set<uint64>& uniqCoors)
{
    vector<exp_t> data;
    uint32 minx, miny;

    try
    {
        hid_t m_file_id;
        hid_t m_dataspace_id;
        hid_t m_dataset_id;
        hsize_t dims[1];
        hid_t attr;

        m_file_id  = H5Fopen(tissueFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

        char expName[128]={0};
        sprintf(expName, "/geneExp/bin1/expression");
        m_dataset_id = H5Dopen(m_file_id, expName, H5P_DEFAULT);
        m_dataspace_id = H5Dget_space(m_dataset_id);
        H5Sget_simple_extent_dims(m_dataspace_id, dims, NULL);
        auto expLen = dims[0];
#ifdef DEBUG
        cout<<"gef expression lines: "<<expLen<<endl;
#endif

        hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(exp_t));
        H5Tinsert(memtype, "x", HOFFSET(exp_t, exp_t::x), H5T_NATIVE_UINT);
        H5Tinsert(memtype, "y", HOFFSET(exp_t, exp_t::y), H5T_NATIVE_UINT);
        H5Tinsert(memtype, "count", HOFFSET(exp_t, exp_t::cnt), H5T_NATIVE_UINT);

        data.resize(expLen);
        H5Dread(m_dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);

        attr = H5Aopen(m_dataset_id, "minX", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_UINT, &minx);
        attr = H5Aopen(m_dataset_id, "minY", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_UINT, &miny);
#ifdef DEBUG
        cout<<"gef offset: "<<minx<<" "<<miny<<endl;
#endif

        H5Aclose(attr);
        H5Tclose(memtype);
        H5Dclose(m_dataset_id);
        H5Sclose(m_dataspace_id);
        H5Fclose(m_file_id);
    }
    catch (...)
    {
        save_error_code(3, "invalid tissue file: " + tissueFile);
        exit(-1);
    }

    // Extract unique coordinates
    for (auto& t : data)
    {
        uint64 coor = ((uint64)(t.x+minx) << 32) + t.y+miny;
        uniqCoors.insert(coor);
    }
}

// Used for hdf5 library compiled with supporting cpp,
// unforunately, we can't use it in cluster
// void parseGefTissueH5CPP(string tissueFile, unordered_set<uint64>& uniqCoors)
// {
//     vector<exp_t> data;
//     uint32 minx, miny;

//     try
//     {
//         H5File* file = new H5File(tissueFile.c_str(), H5F_ACC_RDONLY);
//         DataSet* dataset = new DataSet(file->openDataSet("/geneExp/bin1/expression"));

//         DataSpace dataspace = dataset->getSpace();
//         hsize_t dims[1];
//         dataspace.getSimpleExtentDims(dims, NULL);

//         Attribute attr = dataset->openAttribute("minX");
//         attr.read(PredType::NATIVE_UINT32, &minx);
//         attr = dataset->openAttribute("minY");
//         attr.read(PredType::NATIVE_UINT32, &miny);

//         CompType mtype(sizeof(exp_t));
//         mtype.insertMember("x", HOFFSET(exp_t, x), PredType::NATIVE_UINT32);
//         mtype.insertMember("y", HOFFSET(exp_t, y), PredType::NATIVE_UINT32);
//         mtype.insertMember("count", HOFFSET(exp_t, cnt), PredType::NATIVE_UINT32);

//         data.resize(dims[0]);
//         dataset->read(&data[0], mtype);
        
//         delete dataset;
//         delete file;
//     }
//     catch (...)
//     {
//         cerr<<"SAW-A80001 invalid tissueFile: "<<tissueFile<<endl;
//         throw std::runtime_error("invalid tissueFile!");
//     }

//     // Extract unique coordinates
//     for (auto& t : data)
//     {
//         uint64 coor = ((uint64)(t.x+minx) << 32) + t.y+miny;
//         uniqCoors.insert(coor);
//     }
// }

void parseGemTissue(string tissueFile, unordered_set<uint64>& uniqCoors)
{
    // Parse header for minx, miny
    uint32 minx = 0, miny = 0;
    auto parseHeader = [&](std::string& line){
        if (line[0] != '#') return false;

        if (libx::startswith(line, "#OffsetX"))
        {
            vector<string> temp;
            libx::split(line, '=', temp);
            minx = std::stoi(temp[1]);
        }
        else if (libx::startswith(line, "#OffsetY"))
        {
            vector<string> temp;
            libx::split(line, '=', temp);
            miny = std::stoi(temp[1]);
        }
        return true;
    };
    if (libx::endswith(tissueFile, ".gz"))
        libx::readGZFile(tissueFile, parseHeader, "!", 0);
    else
        libx::readFile(tissueFile, parseHeader, "!", 0);

    // Parse whole data
    vector<uint64> tempCoors;
    auto parseData = [&](std::string& line){
        // Extract unique coordinates=
        string gene;
        int x, y, cnt;
        libx::split(line, '\t', gene, x, y, cnt);
        uint64 coor = ((uint64)(x+minx) << 32) + y+miny;
        tempCoors.emplace_back(coor);

        return true;
    };
    if (libx::endswith(tissueFile, ".gz"))
        libx::readGZFile(tissueFile, parseData, "#", 1);
    else
        libx::readFile(tissueFile, parseData, "#", 1);

    for (auto& coor : tempCoors)
        uniqCoors.insert(coor);

}

bool extractTissueCoor(string tissueFile, unordered_set<uint64>& tissueSampleCoors)
{
    unordered_set<uint64> tissueRawCoors;
    if (libx::endswith(tissueFile, ".gef"))
        parseGefTissue(tissueFile, tissueRawCoors);
    else if (libx::endswith(tissueFile, ".gem.gz") || libx::endswith(tissueFile, ".gem"))
        parseGemTissue(tissueFile, tissueRawCoors);
    
    // Sampleing uniq corrdinates with bin200
    unordered_set<uint64> binCoors;
    for (auto& coor : tissueRawCoors)
    {
        uint32 y = coor & 0xFFFFFFFF;
        uint32 x = (coor >> 32) & 0xFFFFFFFF;
        uint64 binCoor = ((uint64)x/BIN_SIZE << 32) + y/BIN_SIZE;
        binCoors.insert(binCoor);
    }

    float sampleRatio = DATA_RATIO;
    int sampleNum = binCoors.size() * sampleRatio;
    if (sampleNum >= MIN_BIN_THRE)
    {
        vector<uint64> vecBinCoors(binCoors.begin(), binCoors.end());
        std::random_shuffle(vecBinCoors.begin(), vecBinCoors.end());
        binCoors.clear();
        binCoors.insert(vecBinCoors.begin(), vecBinCoors.begin()+sampleNum);
#ifdef DEBUG
        cout<<"gef bin  coordinates number: "<<vecBinCoors.size()<<endl;
#endif
    }

    tissueSampleCoors.clear();
    for (auto& coor : tissueRawCoors)
    {
        uint32 y = coor & 0xFFFFFFFF;
        uint32 x = (coor >> 32) & 0xFFFFFFFF;
        uint64 binCoor = ((uint64)x/BIN_SIZE << 32) + y/BIN_SIZE;
        if (binCoors.count(binCoor) != 0)
            tissueSampleCoors.insert(coor);
    }

#ifdef DEBUG
    cout<<"gef raw  coordinates number: "<<tissueRawCoors.size()<<endl;
    cout<<"gef uniq coordinates number: "<<tissueSampleCoors.size()<<endl;
#endif

    return true;
}

uint64 extractExp(string expFile, unordered_set<uint64>& tissueSampleCoors, vector<ExpLine>& data)
{
    uint64 totalReads = 0;

    // Filter raw data by uniq coordinates after tissuecut
    auto readLine = [&](string& line) {
        uint64 cx, cy, gene, mid, cnt;
        libx::split(line, ' ', cy, cx, gene, mid, cnt);
        uint64 coor = ((uint64)cx << 32) + cy;
        if (tissueSampleCoors.count(coor) != 0)
        {
            ExpLine t(cx, cy, gene, mid);
            data.insert(data.end(), cnt, t);
        }
        totalReads += cnt;
        return true;
    };
    libx::readFile(expFile, readLine, "y");

#ifdef DEBUG
    cout<<"filtered and sampled gene expression lines: "<<data.size()<<endl;
    cout<<"total reads under tissue: "<<totalReads<<endl;
#endif

    // Check is there exists data after filtered by tissue.gef
    if (data.empty())
    {
        save_error_code(4, "no data left after filter by coordinates.");
        exit(-1);
    }

    return totalReads;
}

uint32 calcMedian(vector<uint32>& vec)
{
    if (vec.empty()) return 0;

    const auto midIter = vec.begin() + vec.size() / 2;
    std::nth_element(vec.begin(), midIter, vec.end());
    if (vec.size() % 2 == 0)
    {
        const auto& leftMidIter = std::max_element(vec.begin(),
            midIter);
        return (*leftMidIter + * midIter) / 2;
    }
    else
    {
        return *midIter;
    }
}

bool stat(SatStatMap& statMap, uint64& total, uint64& uniq, uint64& med)
{
    uint64 _total = 0, _uniq = 0;
    vector<uint32> geneNumList;
    unordered_set<uint32> uniqGeneSet;

    // structure: {Coor : {(Gene,MID) : cnt}}
    for (auto& [_, p] : statMap)
    {
        uniqGeneSet.clear();
        for (auto& [t, cnt] : p)
        {
            _total += cnt;
            uniqGeneSet.insert(std::get<0>(t));
        }
        _uniq += p.size();
        geneNumList.push_back(uniqGeneSet.size());
    }

    total = _total;
    uniq = _uniq;
    med = calcMedian(geneNumList);

    return true;
}

bool sample(vector<ExpLine>& data, uint64 totalReads, string satFile)
{
    SatStatMap rawStat;
    SatStatMap binStat;
    MIDStatMap midStat;
    stringstream ss;
    ss << "sample\tbar_x\tbar_y1\tbar_y2\tbar_umi\tbin_x\tbin_y1\tbin_y2\tbin_umi\n";
    
    std::random_shuffle(data.begin(), data.end());
    uint64 pos = 0;
    libx::Timer timer;
    int len = sizeof(SAMPLE_RATIOS) / sizeof(float);
    for (int i = 0; i < len; ++i)
    {
        float ratio = SAMPLE_RATIOS[i];
        ss<<ratio;
        uint64 endPos = data.size() * ratio;
        if (i == len-1)
            endPos = data.size();

        while (pos++ < endPos)
        {
            auto& t = data[pos];
            GeneMID gm(t.gene, t.mid);

            // collect raw data(bin1)
            uint64 coor = ((uint64)t.cx << 32) + t.cy;
            rawStat[coor][gm]++;

            // collect bin data(bin200)
            uint64 binCoor = ((uint64)t.cx/BIN_SIZE << 32) + t.cy/BIN_SIZE;
            binStat[binCoor][gm]++;

            // collect unique mid in bin200
            CoorGeneMID cgm(t.cx, t.cy, t.gene, t.mid);
            midStat[binCoor].insert(cgm);
        }

        uint64 realReads = totalReads * ratio;
        if (i == len-1)
            realReads = totalReads;
        // synchronize
        uint64 total, uniq, geMed;
        float sat;
        stat(rawStat, total, uniq, geMed);
        sat = total > 0 ? 1-uniq*1.0/total : 0;
        uint64 realUniq = total > 0 ? uniq * realReads / total : 0;
        ss<<"\t"<<realReads<<"\t"<<sat<<"\t"<<geMed<<"\t"<<realUniq;

        stat(binStat, total, uniq, geMed);
        sat = total > 0 ? 1-uniq*1.0/total : 0;
        uint64 sumCnt = 0;
        for (auto& [_, uniqSet] : midStat)
            sumCnt += uniqSet.size();
        uint64 aveCnt = midStat.empty() ? 0 : sumCnt / midStat.size();
        ss<<"\t"<<realReads<<"\t"<<sat<<"\t"<<geMed<<"\t"<<aveCnt;
        ss<<endl;

#ifdef DEBUG
        cout<<"ratio: "<<ratio<<" reads: "<<endPos<<" timer(s): "<<timer.toc()<<endl;
#endif
    }

    libx::writeFile(ss.str(), satFile);

    return true;
}

bool saturation(string tissueGefFile, string expFile, string satFile)
{
    libx::Timer timer;
    vector<ExpLine> data;
    uint64 totalReads;

    {
        unordered_set<uint64> tissueSampleCoors;
        extractTissueCoor(tissueGefFile, tissueSampleCoors);

#ifdef DEBUG
        cout<<"time(s) of read tissue.gef: "<<timer.toc()<<endl;
#endif

        totalReads = extractExp(expFile, tissueSampleCoors, data);

#ifdef DEBUG
        cout<<"time(s) of read raw gene exp file: "<<timer.toc()<<endl;
        // cout<<"maximum memory cost: "<<libx::getSelfMemory()<<endl;
#endif
    }

    sample(data, totalReads, satFile);

#ifdef DEBUG
    cout<<"time(s) of calculation: "<<timer.toc()<<endl;
#endif

    return true;
}

bool plot(string satFile, float ratio, string outDir)
{
    libx::Timer timer;
    fs::path scriptPath = binPath / "plot.pyc";
    std::string cmd("python3 -W ignore " + scriptPath.string() + " " + satFile +
                     " " + to_string(ratio) + " " + outDir + " 2>&1");
#ifdef DEBUG
    cout<<"run plot command: "<<cmd<<endl;
#endif

    std::string out;
    int rtn = libx::subprocess(cmd, out);
    if (rtn != 0)
    {
        cerr<<out<<endl;
        return false;
    }

#ifdef DEBUG
    cout<<"time(s) of plot: "<<timer.toc()<<endl;
#endif

    return true;
}

float prepare(vector<string>& annoFiles, vector<string>& mapFiles)
{
    // Extract the thrid value(annotated reads number) of thrid line in annotation stat file
    uint64 annoRead = 0;
    auto parseAnnoFile = [&](string& line) {
        vector<string> temp;
        libx::split(line, '\t', temp);
        annoRead = std::stoul(temp[2]);
        return false;
    };
    uint64 totalAnnoRead = 0;
    for (auto& annoFile : annoFiles)
    {
        annoRead = 0;
        libx::readFile(annoFile, parseAnnoFile, "#", 1);
        totalAnnoRead += annoRead;
    }

    // Extract the total reads in bcStar stat file
    uint64 mapRead = 0;
    auto parseMapFile = [&](string& line) {
        if (libx::startswith(line, "total_reads"))
        {
            vector<string> temp;
            libx::split(line, '\t', temp);
            mapRead = std::stoul(temp[1]);
            return false;
        }
        else
        {
            return true;
        }
    };
    uint64 totalMapRead = 0;
    for (auto& mapFile : mapFiles)
    {
        mapRead = 0;
        libx::readFile(mapFile, parseMapFile);
        totalMapRead += mapRead;
    }

    if (totalMapRead == 0)
    {
        save_error_code(5, "total map reads is 0, please check file format from --bcstat");
        exit(-1);
    }
    else if (totalMapRead < totalAnnoRead)
    {
        save_error_code(6, "map reads less than annotated reads.");
        exit(-1);
    }

    float res = totalAnnoRead * 1.0 / totalMapRead;
    
#ifdef DEBUG
    cout<<"anno reads: "<<totalAnnoRead<<" map reads: "<<totalMapRead<<" ratio: "<<res<<endl;
#endif

    return res;
}

int main(int argc, char** argv)
{
    // Parse arguments
    if (!parseArgs(argc, argv))
    {
        // cerr<<"Fail to parse arguments!"<<endl;
        return -1;
    }
    
    // Prepare
    float ratio = prepare(arguments.annoStatFiles, arguments.bcMapStatFiles);

    // Generate saturation file
    if (!saturation(arguments.tissueFile, arguments.expFile, arguments.outSatFile))
    {
        save_error_code(8, "fail to generate saturation file.");
        exit(-1);
    }

    // Plot
    if (!plot(arguments.outSatFile, ratio, arguments.outDir))
    {
        save_error_code(7, "plot error.");
        exit(-1);
    }

    return 0;
}
