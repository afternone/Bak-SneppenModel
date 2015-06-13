#ifndef AVALANCHE_H
#define AVALANCHE_H

#include "heap1.h"
#include <cfloat> // DBL_MAX
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <functional>
#include <limits>
#include <cmath>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
using namespace std;

typedef pair<int,int> ic;
struct Distance
{
    double r;
    int count;
};
typedef pair<double,int> rbar;
//const double minr = 0.0;
//const double maxr = 1.0;

boost::mt19937 rng1(43);
boost::mt19937 rng2(44);
static boost::uniform_01<boost::mt19937> ui(rng1);
//static boost::uniform_01<boost::mt19937> dis(rng2);
boost::exponential_distribution<>expdis;
//boost::variate_generator<boost::mt19937 &, boost::uniform_01<> > ui(gen, ui01);
boost::variate_generator<boost::mt19937 &, boost::exponential_distribution<> > dis(rng2, expdis);

template <typename T>
string num2str(const T &i) {
  string s;
  stringstream ss(s);
  ss << i;
  return ss.str();
}

/*************** 1d *******************/
void ava1d(double f0, double alpha, int sco, int N)
{
    MinHeap<double> h; // the heap to store nodes
    vector<int> hk;
    boost::unordered_map<int, rbar> sr;
    boost::unordered_map<int, rbar>::iterator sr_it;
    boost::unordered_map<int, ic> kh;
    boost::unordered_map<int, ic>::iterator kh_it;
	double r,meanx=0,rds_gyr,sumx; // a temporary random number

	string ofname=string("1_")+num2str(alpha)+string("_")+num2str(f0)+string("_s_ncov_R");
	ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        hk.clear();
		kh.clear();
		h.clear();
		int s=0;
		int tph=0;
		int min_node=0,nei_node;

        // set the active node
        hk.push_back(min_node);
		kh.insert(make_pair(min_node,make_pair(h.push(0.0),0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min node
            sr_it=sr.find(s);
			if(sr_it!=sr.end()){
                (sr_it->second).first+=abs(min_node);
                ++(sr_it->second).second;
			}
			else{
                sr.insert(make_pair(s,make_pair(abs(min_node),1)));
			}
			++kh[min_node].second;
            // ***********update neighbors*************//
			// neighbor 1
			if(ui()<=alpha){
				nei_node=min_node-1;
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 2
			if(ui()<=alpha){
				nei_node=min_node+1;
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			// reset
            tph=h.toph();
            min_node=hk[tph];
		}

		sumx=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            sumx+=(kh_it->first) * (kh_it->second).second;
		}
		meanx=sumx/s/3.0;
		rds_gyr=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            rds_gyr+=((kh_it->first)-meanx)*((kh_it->first)-meanx);
		}
		fout<<s<<"\t"<<hk.size()<<"\t"<<sqrt(rds_gyr/hk.size())<<endl;
	}
	fout.close();
	ofname=string("1_")+num2str(alpha)+string("_")+num2str(f0)+string("_sr");
	fout.open(ofname.c_str());
	for(sr_it=sr.begin();sr_it!=sr.end();++sr_it){
        fout<<sr_it->first<<"\t"<<(sr_it->second).first/(sr_it->second).second<<endl;
	}
	fout.close();
}

/*************** 2d *******************/
struct point2d
{
    int x;
    int y;
public:
    point2d() : x(0), y(0) {}
    point2d(int _x, int _y) : x(_x), y(_y) {}

    bool operator==(point2d const& other) const
    {
        return x == other.x && y == other.y;
    }
    friend std::size_t hash_value(point2d const& p)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.x);
        boost::hash_combine(seed, p.y);

        return seed;
    }
};

struct site2d
{
    int x;
    int y;
public:
    site2d() : x(0), y(0) {}
    site2d(int _x, int _y) : x(_x), y(_y) {}

    bool operator==(site2d const& other) const
    {
        return x == other.x && y == other.y;
    }
    friend std::size_t hash_value(site2d const& p)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.x);
        boost::hash_combine(seed, p.y);

        return seed;
    }
};

void bs2d(double f0, double alpha, int sco, int N)
{
    MinHeap<double> h; // the heap to store sites
    vector<site2d> sites; // from handle to site, sites[handle]
    vector<int> weight; // weight for each sites, i.e. covered times
    Distance *rs=new Distance[sco]; // to store positive distance r relative to the origin at update step s
    boost::unordered_map<site2d, size_t> handles; // a map to from node to handle
    boost::unordered_map<site2d, size_t>::iterator handles_it;
	double r; // a temporary random number
	double sumx,sumy, sumweight; // temporary variable for calculate center of mass
	double meanx=0,meany=0; // center of mass of a single avalanche
	double rds_gyr; // radius of gyration

	Distance tempr;
	tempr.r=0.0;
	tempr.count=0;
	for(int i=0; i<sco; ++i)
	{
	    rs[i]=tempr;
	}

	string ofname=string("2_")+num2str(alpha)+string("_")+num2str(f0)+string("_snR");
	ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        sites.clear();
		handles.clear();
		h.clear();
		int s=0;
		int tph=0;
		site2d min_site(0,0), nei_site;

        // set the active node
        sites.push_back(min_site);
        weight.push_back(1);
		handles.insert(make_pair(min_site, h.push(0.0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min node

			// record the instances of having activity at a positive distance r relative to the origin at update step s
			rs[s].r+=abs(min_site.x)+abs(min_site.y);
			rs[s].count++;

            // ***********update neighbors*************//
			// neighbor 1
			if(ui()<=alpha){
                nei_site=site2d(min_site.x-1,min_site.y);
				handles_it=handles.find(nei_site);
				r=dis();
				if(handles_it!=handles.end()){
					h.update(handles_it->second,r);
					weight[handles_it->second]++;
				}
				else{
                    sites.push_back(nei_site);
                    weight.push_back(1);
                    handles.insert(make_pair(nei_site, h.push(r)));
				}
			}

			//neighbor 2
			if(ui()<=alpha){
                nei_site=site2d(min_site.x+1,min_site.y);
				handles_it=handles.find(nei_site);
				r=dis();
				if(handles_it!=handles.end()){
					h.update(handles_it->second,r);
					weight[handles_it->second]++;
				}
				else{
                    sites.push_back(nei_site);
                    weight.push_back(1);
                    handles.insert(make_pair(nei_site, h.push(r)));
				}
			}

			//neighbor 3
			if(ui()<=alpha){
                nei_site=site2d(min_site.x,min_site.y-1);
				handles_it=handles.find(nei_site);
				r=dis();
				if(handles_it!=handles.end()){
					h.update(handles_it->second,r);
					weight[handles_it->second]++;
				}
				else{
                    sites.push_back(nei_site);
                    weight.push_back(1);
                    handles.insert(make_pair(nei_site, h.push(r)));
				}
			}

			//neighbor 4
			if(ui()<=alpha){
                nei_site=site2d(min_site.x,min_site.y+1);
				handles_it=handles.find(nei_site);
				r=dis();
				if(handles_it!=handles.end()){
					h.update(handles_it->second,r);
					weight[handles_it->second]++;
				}
				else{
                    sites.push_back(nei_site);
                    weight.push_back(1);
                    handles.insert(make_pair(nei_site, h.push(r)));
				}
			}

			// reset
            tph=h.toph();
            min_site=sites[tph];
		}

		sumx=0; sumy=0; sumweight=0;
		for(size_t j=0; j<sites.size(); ++j)
        {
            sumx+=sites[j].x*weight[j];
            sumy+=sites[j].y*weight[j];
            sumweight+=weight[j];
        }
		meanx=sumx/sumweight;
		meany=sumy/sumweight;
		rds_gyr=0;
		for(size_t j=0; j<sites.size(); ++j){
            rds_gyr+=weight[j]*((sites[j].x-meanx)*(sites[j].x-meanx)+(sites[j].y-meany)*(sites[j].y-meany));
		}
		fout<<s<<"\t"<<sites.size()<<"\t"<<sqrt(rds_gyr/sumweight)<<endl;
	}
	fout.close();
	ofname=string("2_")+num2str(alpha)+string("_")+num2str(f0)+string("_sr");
	fout.open(ofname.c_str());
	for(size_t j=0; j<sco; ++j)
    {
        if(rs[j].count>0)
        {
            fout<<j<<"\t"<<rs[j].r/rs[j].count<<endl;
        }
    }
	fout.close();
	delete[] rs;
}

void ava2d(double f0, double alpha, int sco, int N)
{
    MinHeap<double> h; // the heap to store nodes
    vector<point2d> hk;
    boost::unordered_map<int, rbar> sr;
    boost::unordered_map<int, rbar>::iterator sr_it;
    boost::unordered_map<point2d, ic> kh;
    boost::unordered_map<point2d, ic>::iterator kh_it;
	double r,meanx=0,meany=0,rds_gyr,sumx,sumy; // a temporary random number

	string ofname=string("2_")+num2str(alpha)+string("_")+num2str(f0)+string("_s_ncov_R");
	ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        hk.clear();
		kh.clear();
		h.clear();
		int s=0;
		int tph=0;
		point2d min_node(0,0), nei_node;

        // set the active node
        hk.push_back(min_node);
		kh.insert(make_pair(min_node,make_pair(h.push(0.0),0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min node
            sr_it=sr.find(s);
			if(sr_it!=sr.end()){
                (sr_it->second).first+=abs(min_node.x)+abs(min_node.y);
                ++(sr_it->second).second;
			}
			else{
                sr.insert(make_pair(s,make_pair(abs(min_node.x)+abs(min_node.y),1)));
			}
			++kh[min_node].second;
            // ***********update neighbors*************//
			// neighbor 1
			if(ui()<=alpha){
                nei_node=point2d(min_node.x-1,min_node.y);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 2
			if(ui()<=alpha){
                nei_node=point2d(min_node.x+1,min_node.y);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 3
			if(ui()<=alpha){
                nei_node=point2d(min_node.x,min_node.y-1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 4
			if(ui()<=alpha){
                nei_node=point2d(min_node.x,min_node.y+1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			// reset
            tph=h.toph();
            min_node=hk[tph];
		}

		sumx=0; sumy=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            sumx+=(kh_it->first).x * (kh_it->second).second;
            sumy+=(kh_it->first).y * (kh_it->second).second;
		}
		meanx=sumx/s/5.0;
		meany=sumy/s/5.0;
		rds_gyr=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            rds_gyr+=((kh_it->first).x-meanx)*((kh_it->first).x-meanx)
            +((kh_it->first).y-meany)*((kh_it->first).y-meany);
		}
		fout<<s<<"\t"<<hk.size()<<"\t"<<sqrt(rds_gyr/hk.size())<<endl;
	}
	fout.close();
	ofname=string("2_")+num2str(alpha)+string("_")+num2str(f0)+string("_sr");
	fout.open(ofname.c_str());
	for(sr_it=sr.begin();sr_it!=sr.end();++sr_it){
        fout<<sr_it->first<<"\t"<<(sr_it->second).first/(sr_it->second).second<<endl;
	}
	fout.close();
}

/*************** 3d *******************/
struct point3d
{
    int x;
    int y;
    int z;
public:
    point3d() : x(0), y(0), z(0) {}
    point3d(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {}

    bool operator==(point3d const& other) const
    {
        return x == other.x && y == other.y && z == other.z;
    }
    friend std::size_t hash_value(point3d const& p)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.x);
        boost::hash_combine(seed, p.y);
        boost::hash_combine(seed, p.z);

        return seed;
    }
};

void ava3d(double f0, double alpha, int sco, int N)
{
    MinHeap<double> h; // the heap to store nodes
    vector<point3d> hk;
    boost::unordered_map<int, rbar> sr;
    boost::unordered_map<int, rbar>::iterator sr_it;
    boost::unordered_map<point3d, ic> kh;
    boost::unordered_map<point3d, ic>::iterator kh_it;
	double r,tempr,meanx=0,meany=0,meanz=0,rds_gyr,sumx,sumy,sumz; // a temporary random number

	string ofname=string("3_")+num2str(alpha)+string("_")+num2str(f0)+string("_s_ncov_R");
	ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        hk.clear();
		kh.clear();
		h.clear();
		int s=0;
		int tph=0;
		point3d min_node, nei_node;

        // set the active node
        hk.push_back(min_node);
		kh.insert(make_pair(min_node,make_pair(h.push(0.0),0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min node
            sr_it=sr.find(s);
            tempr=abs(min_node.x)+abs(min_node.y)+abs(min_node.z);
			if(sr_it!=sr.end()){
                (sr_it->second).first+=tempr;
                ++(sr_it->second).second;
			}
			else{
                sr.insert(make_pair(s,make_pair(tempr,1)));
			}
			++kh[min_node].second;
            // ***********update neighbors*************//
			// neighbor 1
			if(ui()<=alpha){
                nei_node=point3d(min_node.x-1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 2
			if(ui()<=alpha){
                nei_node=point3d(min_node.x+1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 3
			if(ui()<=alpha){
                nei_node=point3d(min_node.x,min_node.y-1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 4
			if(ui()<=alpha){
                nei_node=point3d(min_node.x,min_node.y+1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 5
			if(ui()<=alpha){
                nei_node=point3d(min_node.x,min_node.y,min_node.z-1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 6
			if(ui()<=alpha){
                nei_node=point3d(min_node.x,min_node.y,min_node.z+1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			// reset
            tph=h.toph();
            min_node=hk[tph];
		}

		sumx=0; sumy=0; sumz=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            sumx+=(kh_it->first).x * (kh_it->second).second;
            sumy+=(kh_it->first).y * (kh_it->second).second;
            sumz+=(kh_it->first).z * (kh_it->second).second;
		}
		meanx=sumx/s/7.0;
		meany=sumy/s/7.0;
		meanz=sumz/s/7.0;
		rds_gyr=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            rds_gyr+=((kh_it->first).x-meanx)*((kh_it->first).x-meanx)
            +((kh_it->first).y-meany)*((kh_it->first).y-meany)
            +((kh_it->first).z-meanz)*((kh_it->first).z-meanz);
		}
		fout<<s<<"\t"<<hk.size()<<"\t"<<sqrt(rds_gyr/hk.size())<<endl;
	}
	fout.close();
	ofname=string("3_")+num2str(alpha)+string("_")+num2str(f0)+string("_sr");
	fout.open(ofname.c_str());
	for(sr_it=sr.begin();sr_it!=sr.end();++sr_it){
        fout<<sr_it->first<<"\t"<<(sr_it->second).first/(sr_it->second).second<<endl;
	}
	fout.close();
}

/*************** 4d *******************/
struct point4d
{
    int w;
    int x;
    int y;
    int z;
public:
    point4d() : w(0), x(0), y(0), z(0) {}
    point4d(int _w, int _x, int _y, int _z) : w(_w), x(_x), y(_y), z(_z) {}

    bool operator==(point4d const& other) const
    {
        return w == other.w && x == other.x && y == other.y && z == other.z;
    }
    friend std::size_t hash_value(point4d const& p)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.w);
        boost::hash_combine(seed, p.x);
        boost::hash_combine(seed, p.y);
        boost::hash_combine(seed, p.z);

        return seed;
    }
};

void ava4d(double f0, double alpha, int sco, int N)
{
    MinHeap<double> h; // the heap to store nodes
    vector<point4d> hk;
    boost::unordered_map<int, rbar> sr;
    boost::unordered_map<int, rbar>::iterator sr_it;
    boost::unordered_map<point4d, ic> kh;
    boost::unordered_map<point4d, ic>::iterator kh_it;
	double r,tempr,meanw=0,meanx=0,meany=0,meanz=0,rds_gyr,sumw,sumx,sumy,sumz; // a temporary random number

	string ofname=string("4_")+num2str(alpha)+string("_")+num2str(f0)+string("_s_ncov_R");
	ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        hk.clear();
		kh.clear();
		h.clear();
		int s=0;
		int tph=0;
		point4d min_node, nei_node;

        // set the active node
        hk.push_back(min_node);
		kh.insert(make_pair(min_node,make_pair(h.push(0.0),0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min node
            sr_it=sr.find(s);
            tempr=abs(min_node.w)+abs(min_node.x)+abs(min_node.y)+abs(min_node.z);
			if(sr_it!=sr.end()){
                (sr_it->second).first+=tempr;
                ++(sr_it->second).second;
			}
			else{
                sr.insert(make_pair(s,make_pair(tempr,1)));
			}
			++kh[min_node].second;
            // ***********update neighbors*************//
			// neighbor 1
			if(ui()<=alpha){
                nei_node=point4d(min_node.w-1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 2
			if(ui()<=alpha){
                nei_node=point4d(min_node.w+1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 3
			if(ui()<=alpha){
                nei_node=point4d(min_node.w,min_node.x-1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 4
			if(ui()<=alpha){
                nei_node=point4d(min_node.w,min_node.x+1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 5
			if(ui()<=alpha){
                nei_node=point4d(min_node.w,min_node.x,min_node.y-1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 6
			if(ui()<=alpha){
                nei_node=point4d(min_node.w,min_node.x,min_node.y+1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 7
			if(ui()<=alpha){
                nei_node=point4d(min_node.w,min_node.x,min_node.y,min_node.z-1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 8
			if(ui()<=alpha){
                nei_node=point4d(min_node.w,min_node.x,min_node.y,min_node.z+1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			// reset
            tph=h.toph();
            min_node=hk[tph];
		}

		sumw=0; sumx=0; sumy=0; sumz=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            sumw+=(kh_it->first).w * (kh_it->second).second;
            sumx+=(kh_it->first).x * (kh_it->second).second;
            sumy+=(kh_it->first).y * (kh_it->second).second;
            sumz+=(kh_it->first).z * (kh_it->second).second;
		}
		meanw=sumw/s/9.0;
		meanx=sumx/s/9.0;
		meany=sumy/s/9.0;
		meanz=sumz/s/9.0;
		rds_gyr=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            rds_gyr+=((kh_it->first).w-meanw)*((kh_it->first).w-meanw)
            +((kh_it->first).x-meanx)*((kh_it->first).x-meanx)
            +((kh_it->first).y-meany)*((kh_it->first).y-meany)
            +((kh_it->first).z-meanz)*((kh_it->first).z-meanz);
		}
		fout<<s<<"\t"<<hk.size()<<"\t"<<sqrt(rds_gyr/hk.size())<<endl;
	}
	fout.close();
	ofname=string("4_")+num2str(alpha)+string("_")+num2str(f0)+string("_sr");
	fout.open(ofname.c_str());
	for(sr_it=sr.begin();sr_it!=sr.end();++sr_it){
        fout<<sr_it->first<<"\t"<<(sr_it->second).first/(sr_it->second).second<<endl;
	}
	fout.close();
}

/*************** 5d *******************/
struct point5d
{
    int v;
    int w;
    int x;
    int y;
    int z;
public:
    point5d() : v(0), w(0), x(0), y(0), z(0) {}
    point5d(int _v, int _w, int _x, int _y, int _z) : v(_v), w(_w), x(_x), y(_y), z(_z) {}

    bool operator==(point5d const& other) const
    {
        return v == other.v && w == other.w && x == other.x && y == other.y && z == other.z;
    }
    friend std::size_t hash_value(point5d const& p)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.v);
        boost::hash_combine(seed, p.w);
        boost::hash_combine(seed, p.x);
        boost::hash_combine(seed, p.y);
        boost::hash_combine(seed, p.z);

        return seed;
    }
};

void ava5d(double f0, double alpha, int sco, int N)
{
    MinHeap<double> h; // the heap to store nodes
    vector<point5d> hk;
    boost::unordered_map<int, rbar> sr;
    boost::unordered_map<int, rbar>::iterator sr_it;
    boost::unordered_map<point5d, ic> kh;
    boost::unordered_map<point5d, ic>::iterator kh_it;
	double r,tempr,meanv=0,meanw=0,meanx=0,meany=0,meanz=0,rds_gyr,sumv,sumw,sumx,sumy,sumz; // a temporary random number

	string ofname=string("5_")+num2str(alpha)+string("_")+num2str(f0)+string("_s_ncov_R");
	ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        hk.clear();
		kh.clear();
		h.clear();
		int s=0;
		int tph=0;
		point5d min_node, nei_node;

        // set the active node
        hk.push_back(min_node);
		kh.insert(make_pair(min_node,make_pair(h.push(0.0),0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min node
            sr_it=sr.find(s);
            tempr=abs(min_node.v)+abs(min_node.w)+abs(min_node.x)+abs(min_node.y)+abs(min_node.z);
			if(sr_it!=sr.end()){
                (sr_it->second).first+=tempr;
                ++(sr_it->second).second;
			}
			else{
                sr.insert(make_pair(s,make_pair(tempr,1)));
			}
			++kh[min_node].second;
            // ***********update neighbors*************//
			// neighbor 1
			if(ui()<=alpha){
                nei_node=point5d(min_node.v-1,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 2
			if(ui()<=alpha){
                nei_node=point5d(min_node.v+1,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 3
			if(ui()<=alpha){
                nei_node=point5d(min_node.v,min_node.w-1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}


			//neighbor 4
			if(ui()<=alpha){
                nei_node=point5d(min_node.v,min_node.w+1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 5
			if(ui()<=alpha){
                nei_node=point5d(min_node.v,min_node.w,min_node.x-1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 6
			if(ui()<=alpha){
                nei_node=point5d(min_node.v,min_node.w,min_node.x+1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 7
			if(ui()<=alpha){
                nei_node=point5d(min_node.v,min_node.w,min_node.x,min_node.y-1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 8
			if(ui()<=alpha){
                nei_node=point5d(min_node.v,min_node.w,min_node.x,min_node.y+1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 9
			if(ui()<=alpha){
                nei_node=point5d(min_node.v,min_node.w,min_node.x,min_node.y,min_node.z-1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 10
			if(ui()<=alpha){
                nei_node=point5d(min_node.v,min_node.w,min_node.x,min_node.y,min_node.z+1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			// reset
            tph=h.toph();
            min_node=hk[tph];
		}

		sumv=0; sumw=0; sumx=0; sumy=0; sumz=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            sumv+=(kh_it->first).v * (kh_it->second).second;
            sumw+=(kh_it->first).w * (kh_it->second).second;
            sumx+=(kh_it->first).x * (kh_it->second).second;
            sumy+=(kh_it->first).y * (kh_it->second).second;
            sumz+=(kh_it->first).z * (kh_it->second).second;
		}
		meanv=sumv/s/11.0;
		meanw=sumw/s/11.0;
		meanx=sumx/s/11.0;
		meany=sumy/s/11.0;
		meanz=sumz/s/11.0;
		rds_gyr=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            rds_gyr+=((kh_it->first).v-meanv)*((kh_it->first).v-meanv)
            +((kh_it->first).w-meanw)*((kh_it->first).w-meanw)
            +((kh_it->first).x-meanx)*((kh_it->first).x-meanx)
            +((kh_it->first).y-meany)*((kh_it->first).y-meany)
            +((kh_it->first).z-meanz)*((kh_it->first).z-meanz);
		}
		fout<<s<<"\t"<<hk.size()<<"\t"<<sqrt(rds_gyr/hk.size())<<endl;
	}
	fout.close();
	ofname=string("5_")+num2str(alpha)+string("_")+num2str(f0)+string("_sr");
	fout.open(ofname.c_str());
	for(sr_it=sr.begin();sr_it!=sr.end();++sr_it){
        fout<<sr_it->first<<"\t"<<(sr_it->second).first/(sr_it->second).second<<endl;
	}
	fout.close();
}


/*************** 6d *******************/
struct point6d
{
    int u;
    int v;
    int w;
    int x;
    int y;
    int z;
public:
    point6d() : u(0), v(0), w(0), x(0), y(0), z(0) {}
    point6d(int _u, int _v, int _w, int _x, int _y, int _z) : u(_u), v(_v), w(_w), x(_x), y(_y), z(_z) {}

    bool operator==(point6d const& other) const
    {
        return u == other.u && v == other.v && w == other.w && x == other.x && y == other.y && z == other.z;
    }
    friend std::size_t hash_value(point6d const& p)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.u);
        boost::hash_combine(seed, p.v);
        boost::hash_combine(seed, p.w);
        boost::hash_combine(seed, p.x);
        boost::hash_combine(seed, p.y);
        boost::hash_combine(seed, p.z);

        return seed;
    }
};

void ava6d(double f0, double alpha, int sco, int N)
{
    MinHeap<double> h; // the heap to store nodes
    vector<point6d> hk;
    boost::unordered_map<int, rbar> sr;
    boost::unordered_map<int, rbar>::iterator sr_it;
    boost::unordered_map<point6d, ic> kh;
    boost::unordered_map<point6d, ic>::iterator kh_it;
	double r,tempr,meanu=0,meanv=0,meanw=0,meanx=0,meany=0,meanz=0,rds_gyr,sumu,sumv,sumw,sumx,sumy,sumz; // a temporary random number

	string ofname=string("6_")+num2str(alpha)+string("_")+num2str(f0)+string("_s_ncov_R");
	ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        hk.clear();
		kh.clear();
		h.clear();
		int s=0;
		int tph=0;
		point6d min_node, nei_node;

        // set the active node
        hk.push_back(min_node);
		kh.insert(make_pair(min_node,make_pair(h.push(0.0),0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min node
            sr_it=sr.find(s);
            tempr=abs(min_node.u)+abs(min_node.v)+abs(min_node.w)+abs(min_node.x)+abs(min_node.y)+abs(min_node.z);
			if(sr_it!=sr.end()){
                (sr_it->second).first+=tempr;
                ++(sr_it->second).second;
			}
			else{
                sr.insert(make_pair(s,make_pair(tempr,1)));
			}
			++kh[min_node].second;
            // ***********update neighbors*************//
			// neighbor 1
			if(ui()<=alpha){
                nei_node=point6d(min_node.u-1,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 2
			if(ui()<=alpha){
                nei_node=point6d(min_node.u+1,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 3
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v-1,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}


			//neighbor 4
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v+1,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 5
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v,min_node.w-1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 6
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v,min_node.w+1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 7
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v,min_node.w,min_node.x-1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 8
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v,min_node.w,min_node.x+1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 9
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v,min_node.w,min_node.x,min_node.y-1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 10
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v,min_node.w,min_node.x,min_node.y+1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 11
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z-1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 12
			if(ui()<=alpha){
                nei_node=point6d(min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z+1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			// reset
            tph=h.toph();
            min_node=hk[tph];
		}

		sumu=0; sumv=0; sumw=0; sumx=0; sumy=0; sumz=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            sumu+=(kh_it->first).u * (kh_it->second).second;
            sumv+=(kh_it->first).v * (kh_it->second).second;
            sumw+=(kh_it->first).w * (kh_it->second).second;
            sumx+=(kh_it->first).x * (kh_it->second).second;
            sumy+=(kh_it->first).y * (kh_it->second).second;
            sumz+=(kh_it->first).z * (kh_it->second).second;
		}
		meanu=sumu/13.0/s;
		meanv=sumv/s/13.0;
		meanw=sumw/s/13.0;
		meanx=sumx/s/13.0;
		meany=sumy/s/13.0;
		meanz=sumz/s/13.0;
		rds_gyr=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            rds_gyr+=((kh_it->first).u-meanu)*((kh_it->first).u-meanu)
            +((kh_it->first).v-meanv)*((kh_it->first).v-meanv)
            +((kh_it->first).w-meanw)*((kh_it->first).w-meanw)
            +((kh_it->first).x-meanx)*((kh_it->first).x-meanx)
            +((kh_it->first).y-meany)*((kh_it->first).y-meany)
            +((kh_it->first).z-meanz)*((kh_it->first).z-meanz);
		}
		fout<<s<<"\t"<<hk.size()<<"\t"<<sqrt(rds_gyr/hk.size())<<endl;
	}
	fout.close();
	ofname=string("6_")+num2str(alpha)+string("_")+num2str(f0)+string("_sr");
	fout.open(ofname.c_str());
	for(sr_it=sr.begin();sr_it!=sr.end();++sr_it){
        fout<<sr_it->first<<"\t"<<(sr_it->second).first/(sr_it->second).second<<endl;
	}
	fout.close();
}

/*************** 7d *******************/
struct point7d
{
    int t;
    int u;
    int v;
    int w;
    int x;
    int y;
    int z;
public:
    point7d() : t(0), u(0), v(0), w(0), x(0), y(0), z(0) {}
    point7d(int _t, int _u, int _v, int _w, int _x, int _y, int _z) : t(_t), u(_u), v(_v), w(_w), x(_x), y(_y), z(_z) {}

    bool operator==(point7d const& other) const
    {
        return t == other.t && u == other.u && v == other.v && w == other.w && x == other.x && y == other.y && z == other.z;
    }
    friend std::size_t hash_value(point7d const& p)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.t);
        boost::hash_combine(seed, p.u);
        boost::hash_combine(seed, p.v);
        boost::hash_combine(seed, p.w);
        boost::hash_combine(seed, p.x);
        boost::hash_combine(seed, p.y);
        boost::hash_combine(seed, p.z);

        return seed;
    }
};

void ava7d(double f0, double alpha, int sco, int N)
{
    MinHeap<double> h; // the heap to store nodes
    vector<point7d> hk;
    boost::unordered_map<int, rbar> sr;
    boost::unordered_map<int, rbar>::iterator sr_it;
    boost::unordered_map<point7d, ic> kh;
    boost::unordered_map<point7d, ic>::iterator kh_it;
	double r,tempr,meant=0,meanu=0,meanv=0,meanw=0,meanx=0,meany=0,meanz=0,rds_gyr,sumt,sumu,sumv,sumw,sumx,sumy,sumz; // a temporary random number

	string ofname=string("7_")+num2str(alpha)+string("_")+num2str(f0)+string("_s_ncov_R");
	ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        hk.clear();
		kh.clear();
		h.clear();
		int s=0;
		int tph=0;
		point7d min_node, nei_node;

        // set the active node
        hk.push_back(min_node);
		kh.insert(make_pair(min_node,make_pair(h.push(0.0),0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min node
            sr_it=sr.find(s);
            tempr=abs(min_node.t)+abs(min_node.u)+abs(min_node.v)+abs(min_node.w)+abs(min_node.x)+abs(min_node.y)+abs(min_node.z);
			if(sr_it!=sr.end()){
                (sr_it->second).first+=tempr;
                ++(sr_it->second).second;
			}
			else{
                sr.insert(make_pair(s,make_pair(tempr,1)));
			}
			++kh[min_node].second;
            // ***********update neighbors*************//
			// neighbor 1
			if(ui()<=alpha){
                nei_node=point7d(min_node.t-1,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 2
			if(ui()<=alpha){
                nei_node=point7d(min_node.t+1,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 3
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u-1,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}


			//neighbor 4
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u+1,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 5
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v-1,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 6
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v+1,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 7
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v,min_node.w-1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 8
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v,min_node.w+1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 9
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v,min_node.w,min_node.x-1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 10
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v,min_node.w,min_node.x+1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 11
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y-1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 12
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y+1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 13
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z-1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 14
			if(ui()<=alpha){
                nei_node=point7d(min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z+1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			// reset
            tph=h.toph();
            min_node=hk[tph];
		}

		sumt=0; sumu=0; sumv=0; sumw=0; sumx=0; sumy=0; sumz=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            sumt+=(kh_it->first).t * (kh_it->second).second;
            sumu+=(kh_it->first).u * (kh_it->second).second;
            sumv+=(kh_it->first).v * (kh_it->second).second;
            sumw+=(kh_it->first).w * (kh_it->second).second;
            sumx+=(kh_it->first).x * (kh_it->second).second;
            sumy+=(kh_it->first).y * (kh_it->second).second;
            sumz+=(kh_it->first).z * (kh_it->second).second;
		}
		meant=sumt/s/15.0;
		meanu=sumu/15.0/s;
		meanv=sumv/s/15.0;
		meanw=sumw/s/15.0;
		meanx=sumx/s/15.0;
		meany=sumy/s/15.0;
		meanz=sumz/s/15.0;
		rds_gyr=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            rds_gyr+=((kh_it->first).t-meant)*((kh_it->first).t-meant)
            +((kh_it->first).u-meanu)*((kh_it->first).u-meanu)
            +((kh_it->first).v-meanv)*((kh_it->first).v-meanv)
            +((kh_it->first).w-meanw)*((kh_it->first).w-meanw)
            +((kh_it->first).x-meanx)*((kh_it->first).x-meanx)
            +((kh_it->first).y-meany)*((kh_it->first).y-meany)
            +((kh_it->first).z-meanz)*((kh_it->first).z-meanz);
		}
		fout<<s<<"\t"<<hk.size()<<"\t"<<sqrt(rds_gyr/hk.size())<<endl;
	}
	fout.close();
	ofname=string("7_")+num2str(alpha)+string("_")+num2str(f0)+string("_sr");
	fout.open(ofname.c_str());
	for(sr_it=sr.begin();sr_it!=sr.end();++sr_it){
        fout<<sr_it->first<<"\t"<<(sr_it->second).first/(sr_it->second).second<<endl;
	}
	fout.close();
}

/*************** 8d *******************/
struct point8d
{
    int s;
    int t;
    int u;
    int v;
    int w;
    int x;
    int y;
    int z;
public:
    point8d() : s(0), t(0), u(0), v(0), w(0), x(0), y(0), z(0) {}
    point8d(int _s, int _t, int _u, int _v, int _w, int _x, int _y, int _z) : s(_s), t(_t), u(_u), v(_v), w(_w), x(_x), y(_y), z(_z) {}

    bool operator==(point8d const& other) const
    {
        return s == other.s && t == other.t && u == other.u && v == other.v && w == other.w && x == other.x && y == other.y && z == other.z;
    }
    friend std::size_t hash_value(point8d const& p)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.s);
        boost::hash_combine(seed, p.t);
        boost::hash_combine(seed, p.u);
        boost::hash_combine(seed, p.v);
        boost::hash_combine(seed, p.w);
        boost::hash_combine(seed, p.x);
        boost::hash_combine(seed, p.y);
        boost::hash_combine(seed, p.z);

        return seed;
    }
};

void ava8d(double f0, double alpha, int sco, int N)
{
    MinHeap<double> h; // the heap to store nodes
    vector<point8d> hk;
    boost::unordered_map<int, rbar> sr;
    boost::unordered_map<int, rbar>::iterator sr_it;
    boost::unordered_map<point8d, ic> kh;
    boost::unordered_map<point8d, ic>::iterator kh_it;
	double r,tempr,means=0,meant=0,meanu=0,meanv=0,meanw=0,meanx=0,meany=0,meanz=0,rds_gyr,sums,sumt,sumu,sumv,sumw,sumx,sumy,sumz; // a temporary random number

	string ofname=string("8_")+num2str(alpha)+string("_")+num2str(f0)+string("_s_ncov_R");
	ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        hk.clear();
		kh.clear();
		h.clear();
		int s=0;
		int tph=0;
		point8d min_node, nei_node;

        // set the active node
        hk.push_back(min_node);
		kh.insert(make_pair(min_node,make_pair(h.push(0.0),0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min node
            sr_it=sr.find(s);
            tempr=abs(min_node.s)+abs(min_node.t)+abs(min_node.u)+abs(min_node.v)+abs(min_node.w)+abs(min_node.x)+abs(min_node.y)+abs(min_node.z);
			if(sr_it!=sr.end()){
                (sr_it->second).first+=tempr;
                ++(sr_it->second).second;
			}
			else{
                sr.insert(make_pair(s,make_pair(tempr,1)));
			}
			++kh[min_node].second;
            // ***********update neighbors*************//
			// neighbor 1
			if(ui()<=alpha){
                nei_node=point8d(min_node.s-1,min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 2
			if(ui()<=alpha){
                nei_node=point8d(min_node.s+1,min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 3
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t-1,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}


			//neighbor 4
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t+1,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 5
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u-1,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 6
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u+1,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 7
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v-1,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 8
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v+1,min_node.w,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 9
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v,min_node.w-1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 10
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v,min_node.w+1,min_node.x,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 11
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v,min_node.w,min_node.x-1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 12
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v,min_node.w,min_node.x+1,min_node.y,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 13
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y-1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 14
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y+1,min_node.z);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 15
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z-1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			//neighbor 16
			if(ui()<=alpha){
                nei_node=point8d(min_node.s,min_node.t,min_node.u,min_node.v,min_node.w,min_node.x,min_node.y,min_node.z+1);
				kh_it=kh.find(nei_node);
				r=dis();
				if(kh_it!=kh.end()){
					h.update((kh_it->second).first,r);
				}
				else{
                    hk.push_back(nei_node);
                    kh.insert(make_pair(nei_node,make_pair(h.push(r),0)));
				}
				++kh[nei_node].second;
			}

			// reset
            tph=h.toph();
            min_node=hk[tph];
		}

		sums=0; sumt=0; sumu=0; sumv=0; sumw=0; sumx=0; sumy=0; sumz=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
			sums+=(kh_it->first).s * (kh_it->second).second;
            sumt+=(kh_it->first).t * (kh_it->second).second;
            sumu+=(kh_it->first).u * (kh_it->second).second;
            sumv+=(kh_it->first).v * (kh_it->second).second;
            sumw+=(kh_it->first).w * (kh_it->second).second;
            sumx+=(kh_it->first).x * (kh_it->second).second;
            sumy+=(kh_it->first).y * (kh_it->second).second;
            sumz+=(kh_it->first).z * (kh_it->second).second;
		}
		means=sums/s/17.0;
		meant=sumt/s/17.0;
		meanu=sumu/17.0/s;
		meanv=sumv/s/17.0;
		meanw=sumw/s/17.0;
		meanx=sumx/s/17.0;
		meany=sumy/s/17.0;
		meanz=sumz/s/17.0;
		rds_gyr=0;
		for(kh_it=kh.begin();kh_it!=kh.end();++kh_it){
            rds_gyr+=((kh_it->first).s-means)*((kh_it->first).s-means)
            +((kh_it->first).t-meant)*((kh_it->first).t-meant)
            +((kh_it->first).u-meanu)*((kh_it->first).u-meanu)
            +((kh_it->first).v-meanv)*((kh_it->first).v-meanv)
            +((kh_it->first).w-meanw)*((kh_it->first).w-meanw)
            +((kh_it->first).x-meanx)*((kh_it->first).x-meanx)
            +((kh_it->first).y-meany)*((kh_it->first).y-meany)
            +((kh_it->first).z-meanz)*((kh_it->first).z-meanz);
		}
		fout<<s<<"\t"<<hk.size()<<"\t"<<sqrt(rds_gyr/hk.size())<<endl;
	}
	fout.close();
	ofname=string("8_")+num2str(alpha)+string("_")+num2str(f0)+string("_sr");
	fout.open(ofname.c_str());
	for(sr_it=sr.begin();sr_it!=sr.end();++sr_it){
        fout<<sr_it->first<<"\t"<<(sr_it->second).first/(sr_it->second).second<<endl;
	}
	fout.close();
}

#endif // AVALANCHE_H