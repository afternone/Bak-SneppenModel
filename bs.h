#ifndef BS_H
#define BS_H

#include "heap.h"
#include <fstream>
#include <cmath>
#include <boost/unordered_map.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/functional/hash.hpp>

// random numbers
typedef boost::mt19937 gen_type;
gen_type gen(43);
static boost::uniform_01<gen_type> ui(gen); // uniform distribution between [0,1)
//static boost::uniform_01<gen_type> dis(gen);
typedef boost::exponential_distribution<> expdis_type;
expdis_type expdis(1.0); // exponential distribution has lambda equal to 1.0
boost::variate_generator<gen_type&, expdis_type> dis(gen, expdis);


// helper function to transform number to string
template <typename T>
std::string num2str(const T &i) {
  std::string s;
  std::stringstream ss(s);
  ss << i;
  return ss.str();
}

/********************* one dimensional Bak-Sneppen model ***********************/
// f0 is  auxiliary parameter
// alpha is the interaction strength between nearest site
// nei is the number of nearest neighbors
// sco is the avalanche cutoff
// N is the total number of avalanche to record
void bs1d(double f0, double alpha, int nei, int sco, int N, bool isotropic=true)
{
    MinHeap<double> h;
    std::vector<int> site_vec;
    std::vector<int> weight_vec;
    boost::unordered_map<int, std::size_t> handle_map;
    boost::unordered_map<int, std::size_t>::iterator handle_it;
    double *r_array=new double[sco];
    int *n_array=new int[sco];
    double r,meanx,rds_gyr,sumx,sumweight;
    int s,rxy,tph,min_site,nei_site;
    std::string suffix;

    for(int i=0; i<sco; ++i)
    {
        r_array[i]=0.0;
        n_array[i]=0;
    }

    if(isotropic)
    {
        suffix=std::string("_snrR.dat");
    }
    else
    {
        suffix=std::string("_anti_snrR.dat");
    }
    std::string ofname=std::string("1_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	std::ofstream fout(ofname.c_str());

	for(int i=0;i<N;i++){

        // initial
        site_vec.clear();
		handle_map.clear();
		h.clear();
		s=0;
		tph=0;
		min_site=0;
		rxy=0;

        // set the active site
        site_vec.push_back(min_site);
        weight_vec.push_back(1);
		handle_map.insert(std::make_pair(min_site, h.push(0.0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco )
        {

			s++;
			h.update(tph,dis()); // update min site

			// record the instances of having activity at a positive distance r relative to the origin at update step s
			r_array[s-1]+=abs(min_site);
			n_array[s-1]++;

            // ***********update neighbors*************//
			// update right neighbors
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha)
                {
                    nei_site=min_site+j;
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site)) rxy=abs(nei_site);
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }


			if(isotropic)
            {
                // update left neighbors if isotropic is true
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha)
                    {
                        nei_site=min_site-j;
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end())
                        {
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else
                        {
                            if(rxy<abs(nei_site)) rxy=abs(nei_site);
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }
            }

			// reset
            tph=h.toph();
            min_site=site_vec[tph];
		} // while

		sumx=0; sumweight=0;
		for(size_t j=0; j<site_vec.size(); ++j)
        {
            sumx+=site_vec[j]*weight_vec[j];
            sumweight+=weight_vec[j];
        }
		meanx=sumx/sumweight;
		rds_gyr=0;
		for(size_t j=0; j<site_vec.size(); ++j){
            rds_gyr+=weight_vec[j]*(site_vec[j]-meanx)*(site_vec[j]-meanx);
		}
		fout<<s<<"\t"<<site_vec.size()<<"\t"<<rxy<<"\t"<<sqrt(rds_gyr/sumweight)<<std::endl;
	} // for
	fout.close();
	if(isotropic)
    {
        suffix=std::string("_sr.dat");
    }
    else
    {
        suffix=std::string("_anti_sr.dat");
    }
	ofname=std::string("1_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	fout.open(ofname.c_str());
	for(size_t j=0; j<sco; ++j)
    {
        if(n_array[j]>0)
        {
            fout<<j+1<<"\t"<<r_array[j]/n_array[j]<<std::endl;
        }
    }
	fout.close();
	delete[] r_array;
	delete[] n_array;
}

/********************* two dimensional Bak-Sneppen model ***********************/
class site2d
{
public:
    int x;
    int y;

    site2d():x(0),y(0){}
    site2d(int _x, int _y) : x(_x), y(_y) {}

    bool operator==(site2d const& other) const
    {
        return x == other.x && y == other.y;
    }

    friend std::size_t hash_value(site2d const& site)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, site.x);
        boost::hash_combine(seed, site.y);
        return seed;
    }
};

void bs2d(double f0, double alpha, int nei, int sco, int N, bool isotropic=true)
{
    MinHeap<double> h; // the heap to store sites
    std::vector<site2d> site_vec; // from handle to site, site_vec[handle]
    std::vector<int> weight_vec; // weight for each sites, i.e. covered times
    double *r_array=new double[sco]; // positive distance r relative to the origin
    int *n_array=new int[sco];
    boost::unordered_map<site2d, size_t> handle_map; // a map to from node to handle
    boost::unordered_map<site2d, size_t>::iterator handle_it;
	double r; // a temporary random number
	double sumx,sumy, sumweight; // temporary variable for calculate center of mass
	double meanx=0,meany=0; // center of mass of a single avalanche
	double rds_gyr; // radius of gyration
    int s,rxy,tph; // avalanche size and top handle
    std::string suffix;

	// initial r_array and n_array
	for(int i=0; i<sco; ++i)
	{
	    r_array[i]=0.0;
	    n_array[i]=0;
	}

	// set the output file name
	if(isotropic)
    {
        suffix=std::string("_snrR.dat");
    }
    else
    {
        suffix=std::string("_anti_snrR.dat");
    }
	std::string ofname=std::string("2_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	std::ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        site_vec.clear();
		handle_map.clear();
		h.clear();
		s=0;
		tph=0;
		site2d min_site(0,0), nei_site;
		rxy=0;

        // set the active site
        site_vec.push_back(min_site);
        weight_vec.push_back(1);
		handle_map.insert(std::make_pair(min_site, h.push(0.0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min site

			// record the instances of having activity at a positive distance r relative to the origin at update step s
			r_array[s-1]+=abs(min_site.x)+abs(min_site.y);
			n_array[s-1]++;

            // ***********update neighbors*************//
			// update right neighbors along dimension 1
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site2d(min_site.x+j,min_site.y);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.x)+abs(nei_site.y)) rxy=abs(nei_site.x)+abs(nei_site.y);
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }


			// update right neighbors along dimension 2
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site2d(min_site.x,min_site.y+j);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.x)+abs(nei_site.y)) rxy=abs(nei_site.x)+abs(nei_site.y);
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

			if(isotropic)
            {
                // update left neighbors along dimension 1
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site2d(min_site.x-j,min_site.y);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.x)+abs(nei_site.y)) rxy=abs(nei_site.x)+abs(nei_site.y);
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update right neighbors along dimension 1
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site2d(min_site.x,min_site.y-j);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.x)+abs(nei_site.y)) rxy=abs(nei_site.x)+abs(nei_site.y);
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }
            }

			// reset
            tph=h.toph();
            min_site=site_vec[tph];
		}

		sumx=0; sumy=0; sumweight=0;
		for(size_t j=0; j<site_vec.size(); ++j)
        {
            sumx+=site_vec[j].x*weight_vec[j];
            sumy+=site_vec[j].y*weight_vec[j];
            sumweight+=weight_vec[j];
        }
		meanx=sumx/sumweight;
		meany=sumy/sumweight;
		rds_gyr=0;
		for(size_t j=0; j<site_vec.size(); ++j){
            rds_gyr+=weight_vec[j]*((site_vec[j].x-meanx)*(site_vec[j].x-meanx)+(site_vec[j].y-meany)*(site_vec[j].y-meany));
		}
		fout<<s<<"\t"<<site_vec.size()<<"\t"<<rxy<<"\t"<<sqrt(rds_gyr/sumweight)<<std::endl;
	} // for
	fout.close();

	if(isotropic)
    {
        suffix=std::string("_sr.dat");
    }
    else
    {
        suffix=std::string("_anti_sr.dat");
    }
	ofname=std::string("2_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	fout.open(ofname.c_str());
	for(size_t j=0; j<sco; ++j)
    {
        if(n_array[j]>0)
        {
            fout<<j+1<<"\t"<<r_array[j]/n_array[j]<<std::endl;
        }
    }
	fout.close();
	delete[] r_array;
	delete[] n_array;
}

/********************* three dimensional Bak-Sneppen model ***********************/
class site3d
{
public:
    int x;
    int y;
    int z;

    site3d():x(0),y(0),z(0){}
    site3d(int _x, int _y, int _z) : x(_x), y(_y), z(_z){}

    bool operator==(site3d const& other) const
    {
        return x == other.x && y == other.y && z == other.z;
    }

    friend std::size_t hash_value(site3d const& site)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, site.x);
        boost::hash_combine(seed, site.y);
        boost::hash_combine(seed, site.z);
        return seed;
    }
};

void bs3d(double f0, double alpha, int nei, int sco, int N, bool isotropic=true)
{
    MinHeap<double> h; // the heap to store sites
    std::vector<site3d> site_vec; // from handle to site, site_vec[handle]
    std::vector<int> weight_vec; // weight for each sites, i.e. covered times
    double *r_array=new double[sco]; // positive distance r relative to the origin
    int *n_array=new int[sco];
    boost::unordered_map<site3d, size_t> handle_map; // a map to from node to handle
    boost::unordered_map<site3d, size_t>::iterator handle_it;
	double r; // a temporary random number
	double sumx,sumy,sumz,sumweight; // temporary variable for calculate center of mass
	double meanx,meany,meanz; // center of mass of a single avalanche
	double rds_gyr; // radius of gyration
    int s,rxy,tph; // avalanche size and top handle
    std::string suffix; // suffix of the output file

	// initial r_array and n_array
	for(int i=0; i<sco; ++i)
	{
	    r_array[i]=0.0;
	    n_array[i]=0;
	}

	// set the output file name
	if(isotropic)
    {
        suffix=std::string("_snrR.dat");
    }
    else
    {
        suffix=std::string("_anti_snrR.dat");
    }
	std::string ofname=std::string("3_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	std::ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        site_vec.clear();
		handle_map.clear();
		h.clear();
		s=0;
		tph=0;
		site3d min_site(0,0,0), nei_site;
		rxy=0;

        // set the active site
        site_vec.push_back(min_site);
        weight_vec.push_back(1);
		handle_map.insert(std::make_pair(min_site, h.push(0.0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min site

			// record the instances of having activity at a positive distance r relative to the origin at update step s
			r_array[s-1]+=abs(min_site.x)+abs(min_site.y)+abs(min_site.z);
			n_array[s-1]++;

            // ***********update neighbors*************//
			// update right neighbors along dimension 1
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site3d(min_site.x+j,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }


			// update right neighbors along dimension 2
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site3d(min_site.x,min_site.y+j,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 3
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site3d(min_site.x,min_site.y,min_site.z+j);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

			if(isotropic)
            {
                // update left neighbors along dimension 1
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site3d(min_site.x-j,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update right neighbors along dimension 2
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site3d(min_site.x,min_site.y-j,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update right neighbors along dimension 3
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site3d(min_site.x,min_site.y,min_site.z-j);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }
            }

			// reset
            tph=h.toph();
            min_site=site_vec[tph];
		}

		sumx=0; sumy=0; sumz=0; sumweight=0;
		for(size_t j=0; j<site_vec.size(); ++j)
        {
            sumx+=site_vec[j].x*weight_vec[j];
            sumy+=site_vec[j].y*weight_vec[j];
            sumz+=site_vec[j].z*weight_vec[j];
            sumweight+=weight_vec[j];
        }
		meanx=sumx/sumweight;
		meany=sumy/sumweight;
		meanz=sumz/sumweight;
		rds_gyr=0;
		for(size_t j=0; j<site_vec.size(); ++j){
            rds_gyr+=weight_vec[j]*((site_vec[j].x-meanx)*(site_vec[j].x-meanx)
                                    +(site_vec[j].y-meany)*(site_vec[j].y-meany)
                                    +(site_vec[j].z-meanz)*(site_vec[j].z-meanz));
		}
		fout<<s<<"\t"<<site_vec.size()<<"\t"<<rxy<<"\t"<<sqrt(rds_gyr/sumweight)<<std::endl;
	} // for
	fout.close();

	if(isotropic)
    {
        suffix=std::string("_sr.dat");
    }
    else
    {
        suffix=std::string("_anti_sr.dat");
    }
	ofname=std::string("3_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	fout.open(ofname.c_str());
	for(size_t j=0; j<sco; ++j)
    {
        if(n_array[j]>0)
        {
            fout<<j+1<<"\t"<<r_array[j]/n_array[j]<<std::endl;
        }
    }
	fout.close();
	delete[] r_array;
	delete[] n_array;
}

/********************* four dimensional Bak-Sneppen model ***********************/
class site4d
{
public:
    int w;
    int x;
    int y;
    int z;

    site4d():w(0),x(0),y(0),z(0){}
    site4d(int _w, int _x, int _y, int _z) : w(_w), x(_x), y(_y), z(_z){}

    bool operator==(site4d const& other) const
    {
        return w == other.w && x == other.x && y == other.y && z == other.z;
    }

    friend std::size_t hash_value(site4d const& site)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, site.w);
        boost::hash_combine(seed, site.x);
        boost::hash_combine(seed, site.y);
        boost::hash_combine(seed, site.z);
        return seed;
    }
};

void bs4d(double f0, double alpha, int nei, int sco, int N, bool isotropic=true)
{
    MinHeap<double> h; // the heap to store sites
    std::vector<site4d> site_vec; // from handle to site, site_vec[handle]
    std::vector<int> weight_vec; // weight for each sites, i.e. covered times
    double *r_array=new double[sco]; // positive distance r relative to the origin
    int *n_array=new int[sco];
    boost::unordered_map<site4d, size_t> handle_map; // a map to from node to handle
    boost::unordered_map<site4d, size_t>::iterator handle_it;
	double r; // a temporary random number
	double sumw,sumx,sumy,sumz,sumweight; // temporary variable for calculate center of mass
	double meanw,meanx,meany,meanz; // center of mass of a single avalanche
	double rds_gyr; // radius of gyration
    int s,rxy,tph; // avalanche size and top handle
    std::string suffix; // suffix of the output file

	// initial r_array and n_array
	for(int i=0; i<sco; ++i)
	{
	    r_array[i]=0.0;
	    n_array[i]=0;
	}

	// set the output file name
	if(isotropic)
    {
        suffix=std::string("_snrR.dat");
    }
    else
    {
        suffix=std::string("_anti_snrR.dat");
    }
	std::string ofname=std::string("4_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	std::ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        site_vec.clear();
		handle_map.clear();
		h.clear();
		s=0;
		tph=0;
		site4d min_site(0,0,0,0), nei_site;
		rxy=0;

        // set the active site
        site_vec.push_back(min_site);
        weight_vec.push_back(1);
		handle_map.insert(std::make_pair(min_site, h.push(0.0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min site

			// record the instances of having activity at a positive distance r relative to the origin at update step s
			r_array[s-1]+=abs(min_site.w)+abs(min_site.x)+abs(min_site.y)+abs(min_site.z);
			n_array[s-1]++;

            // ***********update neighbors*************//
			// update right neighbors along dimension 1
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site4d(min_site.w+j,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }


			// update right neighbors along dimension 2
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site4d(min_site.w,min_site.x+j,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 3
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site4d(min_site.w,min_site.x,min_site.y+j,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 4
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site4d(min_site.w,min_site.x,min_site.y,min_site.z+j);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

			if(isotropic)
            {
                // update left neighbors along dimension 1
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site4d(min_site.w-j,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 2
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site4d(min_site.w,min_site.x-j,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 3
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site4d(min_site.w,min_site.x,min_site.y-j,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 4
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site4d(min_site.w,min_site.x,min_site.y,min_site.z-j);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }
            }

			// reset
            tph=h.toph();
            min_site=site_vec[tph];
		}

		sumw=0; sumx=0; sumy=0; sumz=0; sumweight=0;
		for(size_t j=0; j<site_vec.size(); ++j)
        {
            sumw+=site_vec[j].w*weight_vec[j];
            sumx+=site_vec[j].x*weight_vec[j];
            sumy+=site_vec[j].y*weight_vec[j];
            sumz+=site_vec[j].z*weight_vec[j];
            sumweight+=weight_vec[j];
        }
        meanw=sumw/sumweight;
		meanx=sumx/sumweight;
		meany=sumy/sumweight;
		meanz=sumz/sumweight;
		rds_gyr=0;
		for(size_t j=0; j<site_vec.size(); ++j){
            rds_gyr+=weight_vec[j]*(+(site_vec[j].w-meanw)*(site_vec[j].w-meanw)
                                    +(site_vec[j].x-meanx)*(site_vec[j].x-meanx)
                                    +(site_vec[j].y-meany)*(site_vec[j].y-meany)
                                    +(site_vec[j].z-meanz)*(site_vec[j].z-meanz));
		}
		fout<<s<<"\t"<<site_vec.size()<<"\t"<<rxy<<"\t"<<sqrt(rds_gyr/sumweight)<<std::endl;
	}
	fout.close();

	if(isotropic)
    {
        suffix=std::string("_sr.dat");
    }
    else
    {
        suffix=std::string("_anti_sr.dat");
    }
	ofname=std::string("4_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	fout.open(ofname.c_str());
	for(size_t j=0; j<sco; ++j)
    {
        if(n_array[j]>0)
        {
            fout<<j+1<<"\t"<<r_array[j]/n_array[j]<<std::endl;
        }
    }
	fout.close();
	delete[] r_array;
	delete[] n_array;
}

/********************* five dimensional Bak-Sneppen model ***********************/
class site5d
{
public:
    int v;
    int w;
    int x;
    int y;
    int z;

    site5d():v(0),w(0),x(0),y(0),z(0){}
    site5d(int _v, int _w, int _x, int _y, int _z) : v(_v), w(_w), x(_x), y(_y), z(_z){}

    bool operator==(site5d const& other) const
    {
        return v==other.v && w==other.w && x==other.x && y==other.y && z==other.z;
    }

    friend std::size_t hash_value(site5d const& site)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, site.v);
        boost::hash_combine(seed, site.w);
        boost::hash_combine(seed, site.x);
        boost::hash_combine(seed, site.y);
        boost::hash_combine(seed, site.z);
        return seed;
    }
};

void bs5d(double f0, double alpha, int nei, int sco, int N, bool isotropic=true)
{
    MinHeap<double> h; // the heap to store sites
    std::vector<site5d> site_vec; // from handle to site, site_vec[handle]
    std::vector<int> weight_vec; // weight for each sites, i.e. covered times
    double *r_array=new double[sco]; // positive distance r relative to the origin
    int *n_array=new int[sco];
    boost::unordered_map<site5d, size_t> handle_map; // a map to from node to handle
    boost::unordered_map<site5d, size_t>::iterator handle_it;
	double r; // a temporary random number
	double sumv,sumw,sumx,sumy,sumz,sumweight; // temporary variable for calculate center of mass
	double meanv,meanw,meanx,meany,meanz; // center of mass of a single avalanche
	double rds_gyr; // radius of gyration
    int s,rxy,tph; // avalanche size and top handle
    std::string suffix; // suffix of the output file

	// initial r_array and n_array
	for(int i=0; i<sco; ++i)
	{
	    r_array[i]=0.0;
	    n_array[i]=0;
	}

	// set the output file name
	if(isotropic)
    {
        suffix=std::string("_snrR.dat");
    }
    else
    {
        suffix=std::string("_anti_snrR.dat");
    }
	std::string ofname=std::string("5_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	std::ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        site_vec.clear();
		handle_map.clear();
		h.clear();
		s=0;
		tph=0;
		site5d min_site(0,0,0,0,0), nei_site;
		rxy=0; // initial r

        // set the active site
        site_vec.push_back(min_site);
        weight_vec.push_back(1);
		handle_map.insert(std::make_pair(min_site, h.push(0.0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min site

			// record the instances of having activity at a positive distance r relative to the origin at update step s
			r_array[s-1]+=abs(min_site.v)+abs(min_site.w)+abs(min_site.x)+abs(min_site.y)+abs(min_site.z);
			n_array[s-1]++;

            // ***********update neighbors*************//
			// update right neighbors along dimension 1
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site5d(min_site.v+j,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }


			// update right neighbors along dimension 2
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site5d(min_site.v,min_site.w+j,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 3
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site5d(min_site.v,min_site.w,min_site.x+j,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 4
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site5d(min_site.v,min_site.w,min_site.x,min_site.y+j,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 5
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site5d(min_site.v,min_site.w,min_site.x,min_site.y,min_site.z+j);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

			if(isotropic)
            {
                // update left neighbors along dimension 1
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site5d(min_site.v-j,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 2
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site5d(min_site.v,min_site.w-j,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 3
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site5d(min_site.v,min_site.w,min_site.x-j,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 4
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site5d(min_site.v,min_site.w,min_site.x,min_site.y-j,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 5
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site5d(min_site.v,min_site.w,min_site.x,min_site.y,min_site.z-j);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }
            }

			// reset
            tph=h.toph();
            min_site=site_vec[tph];
		}

		sumv=0; sumw=0; sumx=0; sumy=0; sumz=0; sumweight=0;
		for(size_t j=0; j<site_vec.size(); ++j)
        {
            sumv+=site_vec[j].v*weight_vec[j];
            sumw+=site_vec[j].w*weight_vec[j];
            sumx+=site_vec[j].x*weight_vec[j];
            sumy+=site_vec[j].y*weight_vec[j];
            sumz+=site_vec[j].z*weight_vec[j];
            sumweight+=weight_vec[j];
        }
        meanv=sumv/sumweight;
        meanw=sumw/sumweight;
		meanx=sumx/sumweight;
		meany=sumy/sumweight;
		meanz=sumz/sumweight;
		rds_gyr=0;
		for(size_t j=0; j<site_vec.size(); ++j){
            rds_gyr+=weight_vec[j]*((site_vec[j].v-meanv)*(site_vec[j].v-meanv)
                                    +(site_vec[j].w-meanw)*(site_vec[j].w-meanw)
                                    +(site_vec[j].x-meanx)*(site_vec[j].x-meanx)
                                    +(site_vec[j].y-meany)*(site_vec[j].y-meany)
                                    +(site_vec[j].z-meanz)*(site_vec[j].z-meanz));
		}
		fout<<s<<"\t"<<site_vec.size()<<"\t"<<rxy<<"\t"<<sqrt(rds_gyr/sumweight)<<std::endl;
	}
	fout.close();

	if(isotropic)
    {
        suffix=std::string("_sr.dat");
    }
    else
    {
        suffix=std::string("_anti_sr.dat");
    }
	ofname=std::string("5_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	fout.open(ofname.c_str());
	for(size_t j=0; j<sco; ++j)
    {
        if(n_array[j]>0)
        {
            fout<<j+1<<"\t"<<r_array[j]/n_array[j]<<std::endl;
        }
    }
	fout.close();
	delete[] r_array;
	delete[] n_array;
}

/********************* six dimensional Bak-Sneppen model ***********************/
class site6d
{
public:
    int u;
    int v;
    int w;
    int x;
    int y;
    int z;

    site6d():u(0),v(0),w(0),x(0),y(0),z(0){}
    site6d(int _u, int _v, int _w, int _x, int _y, int _z) : u(_u), v(_v), w(_w), x(_x), y(_y), z(_z){}

    bool operator==(site6d const& other) const
    {
        return u==other.u && v==other.v && w==other.w && x==other.x && y==other.y && z==other.z;
    }

    friend std::size_t hash_value(site6d const& site)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, site.u);
        boost::hash_combine(seed, site.v);
        boost::hash_combine(seed, site.w);
        boost::hash_combine(seed, site.x);
        boost::hash_combine(seed, site.y);
        boost::hash_combine(seed, site.z);
        return seed;
    }
};

void bs6d(double f0, double alpha, int nei, int sco, int N, bool isotropic=true)
{
    MinHeap<double> h; // the heap to store sites
    std::vector<site6d> site_vec; // from handle to site, site_vec[handle]
    std::vector<int> weight_vec; // weight for each sites, i.e. covered times
    double *r_array=new double[sco]; // positive distance r relative to the origin
    int *n_array=new int[sco];
    boost::unordered_map<site6d, size_t> handle_map; // a map to from node to handle
    boost::unordered_map<site6d, size_t>::iterator handle_it;
	double r; // a temporary random number
	double sumu,sumv,sumw,sumx,sumy,sumz,sumweight; // temporary variable for calculate center of mass
	double meanu,meanv,meanw,meanx,meany,meanz; // center of mass of a single avalanche
	double rds_gyr; // radius of gyration
    int s,rxy,tph; // avalanche size and top handle
    std::string suffix; // suffix of the output file

	// initial r_array and n_array
	for(int i=0; i<sco; ++i)
	{
	    r_array[i]=0.0;
	    n_array[i]=0;
	}

	// set the output file name
	if(isotropic)
    {
        suffix=std::string("_snrR.dat");
    }
    else
    {
        suffix=std::string("_anti_snrR.dat");
    }
	std::string ofname=std::string("6_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	std::ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        site_vec.clear();
		handle_map.clear();
		h.clear();
		s=0;
		tph=0;
		site6d min_site(0,0,0,0,0,0), nei_site;
		rxy=0;

        // set the active site
        site_vec.push_back(min_site);
        weight_vec.push_back(1);
		handle_map.insert(std::make_pair(min_site, h.push(0.0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min site

			// record the instances of having activity at a positive distance r relative to the origin at update step s
			r_array[s-1]+=abs(min_site.u)+abs(min_site.v)+abs(min_site.w)+abs(min_site.x)+abs(min_site.y)+abs(min_site.z);
			n_array[s-1]++;

            // ***********update neighbors*************//
			// update right neighbors along dimension 1
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site6d(min_site.u+j,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }


			// update right neighbors along dimension 2
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site6d(min_site.u,min_site.v+j,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 3
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site6d(min_site.u,min_site.v,min_site.w+j,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 4
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site6d(min_site.u,min_site.v,min_site.w,min_site.x+j,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 5
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site6d(min_site.u,min_site.v,min_site.w,min_site.x,min_site.y+j,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 6
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site6d(min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z+j);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

			if(isotropic)
            {
                // update left neighbors along dimension 1
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site6d(min_site.u-j,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 2
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site6d(min_site.u,min_site.v-j,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 3
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site6d(min_site.u,min_site.v,min_site.w-j,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 4
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site6d(min_site.u,min_site.v,min_site.w,min_site.x-j,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 5
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site6d(min_site.u,min_site.v,min_site.w,min_site.x,min_site.y-j,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 6
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site6d(min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z-j);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }
            }

			// reset
            tph=h.toph();
            min_site=site_vec[tph];
		}

		sumu=0; sumv=0; sumw=0; sumx=0; sumy=0; sumz=0; sumweight=0;
		for(size_t j=0; j<site_vec.size(); ++j)
        {
            sumu+=site_vec[j].u*weight_vec[j];
            sumv+=site_vec[j].v*weight_vec[j];
            sumw+=site_vec[j].w*weight_vec[j];
            sumx+=site_vec[j].x*weight_vec[j];
            sumy+=site_vec[j].y*weight_vec[j];
            sumz+=site_vec[j].z*weight_vec[j];
            sumweight+=weight_vec[j];
        }
        meanu=sumu/sumweight;
        meanv=sumv/sumweight;
        meanw=sumw/sumweight;
		meanx=sumx/sumweight;
		meany=sumy/sumweight;
		meanz=sumz/sumweight;
		rds_gyr=0;
		for(size_t j=0; j<site_vec.size(); ++j){
            rds_gyr+=weight_vec[j]*((site_vec[j].u-meanu)*(site_vec[j].u-meanu)
                                    +(site_vec[j].v-meanv)*(site_vec[j].v-meanv)
                                    +(site_vec[j].w-meanw)*(site_vec[j].w-meanw)
                                    +(site_vec[j].x-meanx)*(site_vec[j].x-meanx)
                                    +(site_vec[j].y-meany)*(site_vec[j].y-meany)
                                    +(site_vec[j].z-meanz)*(site_vec[j].z-meanz));
		}
		fout<<s<<"\t"<<site_vec.size()<<"\t"<<rxy<<"\t"<<sqrt(rds_gyr/sumweight)<<std::endl;
	}
	fout.close();

	if(isotropic)
    {
        suffix=std::string("_sr.dat");
    }
    else
    {
        suffix=std::string("_anti_sr.dat");
    }
	ofname=std::string("6_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	fout.open(ofname.c_str());
	for(size_t j=0; j<sco; ++j)
    {
        if(n_array[j]>0)
        {
            fout<<j+1<<"\t"<<r_array[j]/n_array[j]<<std::endl;
        }
    }
	fout.close();
	delete[] r_array;
	delete[] n_array;
}

/********************* seven dimensional Bak-Sneppen model ***********************/
class site7d
{
public:
    int t;
    int u;
    int v;
    int w;
    int x;
    int y;
    int z;

    site7d():t(0),u(0),v(0),w(0),x(0),y(0),z(0){}
    site7d(int _t, int _u, int _v, int _w, int _x, int _y, int _z) : t(_t), u(_u), v(_v), w(_w), x(_x), y(_y), z(_z){}

    bool operator==(site7d const& other) const
    {
        return t==other.t && u==other.u && v==other.v && w==other.w && x==other.x && y==other.y && z==other.z;
    }

    friend std::size_t hash_value(site7d const& site)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, site.t);
        boost::hash_combine(seed, site.u);
        boost::hash_combine(seed, site.v);
        boost::hash_combine(seed, site.w);
        boost::hash_combine(seed, site.x);
        boost::hash_combine(seed, site.y);
        boost::hash_combine(seed, site.z);
        return seed;
    }
};

void bs7d(double f0, double alpha, int nei, int sco, int N, bool isotropic=true)
{
    MinHeap<double> h; // the heap to store sites
    std::vector<site7d> site_vec; // from handle to site, site_vec[handle]
    std::vector<int> weight_vec; // weight for each sites, i.e. covered times
    double *r_array=new double[sco]; // positive distance r relative to the origin
    int *n_array=new int[sco];
    boost::unordered_map<site7d, size_t> handle_map; // a map to from node to handle
    boost::unordered_map<site7d, size_t>::iterator handle_it;
	double r; // a temporary random number
	double sumt,sumu,sumv,sumw,sumx,sumy,sumz,sumweight; // temporary variable for calculate center of mass
	double meant,meanu,meanv,meanw,meanx,meany,meanz; // center of mass of a single avalanche
	double rds_gyr; // radius of gyration
    int s,rxy,tph; // avalanche size and top handle
    std::string suffix; // suffix of the output file

	// initial r_array and n_array
	for(int i=0; i<sco; ++i)
	{
	    r_array[i]=0.0;
	    n_array[i]=0;
	}

	// set the output file name
	if(isotropic)
    {
        suffix=std::string("_snrR.dat");
    }
    else
    {
        suffix=std::string("_anti_snrR.dat");
    }
	std::string ofname=std::string("7_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	std::ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        site_vec.clear();
		handle_map.clear();
		h.clear();
		s=0;
		tph=0;
		site7d min_site(0,0,0,0,0,0,0), nei_site;
		rxy=0;

        // set the active site
        site_vec.push_back(min_site);
        weight_vec.push_back(1);
		handle_map.insert(std::make_pair(min_site, h.push(0.0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min site

			// record the instances of having activity at a positive distance r relative to the origin at update step s
			r_array[s-1]+=abs(min_site.t)+abs(min_site.u)+abs(min_site.v)+abs(min_site.w)+abs(min_site.x)+abs(min_site.y)+abs(min_site.z);
			n_array[s-1]++;

            // ***********update neighbors*************//
			// update right neighbors along dimension 1
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site7d(min_site.t+j,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }


			// update right neighbors along dimension 2
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site7d(min_site.t,min_site.u+j,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 3
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site7d(min_site.t,min_site.u,min_site.v+j,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 4
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site7d(min_site.t,min_site.u,min_site.v,min_site.w+j,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 5
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site7d(min_site.t,min_site.u,min_site.v,min_site.w,min_site.x+j,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 6
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site7d(min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y+j,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 7
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site7d(min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z+j);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

			if(isotropic)
            {
                // update left neighbors along dimension 1
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site7d(min_site.t-j,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 2
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site7d(min_site.t,min_site.u-j,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 3
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site7d(min_site.t,min_site.u,min_site.v-j,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 4
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site7d(min_site.t,min_site.u,min_site.v,min_site.w-j,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 5
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site7d(min_site.t,min_site.u,min_site.v,min_site.w,min_site.x-j,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 6
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site7d(min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y-j,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 7
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site7d(min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z-j);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }
            }

			// reset
            tph=h.toph();
            min_site=site_vec[tph];
		}

		sumt=0; sumu=0; sumv=0; sumw=0; sumx=0; sumy=0; sumz=0; sumweight=0;
		for(size_t j=0; j<site_vec.size(); ++j)
        {
            sumt+=site_vec[j].t*weight_vec[j];
            sumu+=site_vec[j].u*weight_vec[j];
            sumv+=site_vec[j].v*weight_vec[j];
            sumw+=site_vec[j].w*weight_vec[j];
            sumx+=site_vec[j].x*weight_vec[j];
            sumy+=site_vec[j].y*weight_vec[j];
            sumz+=site_vec[j].z*weight_vec[j];
            sumweight+=weight_vec[j];
        }
        meant=sumt/sumweight;
        meanu=sumu/sumweight;
        meanv=sumv/sumweight;
        meanw=sumw/sumweight;
		meanx=sumx/sumweight;
		meany=sumy/sumweight;
		meanz=sumz/sumweight;
		rds_gyr=0;
		for(size_t j=0; j<site_vec.size(); ++j){
            rds_gyr+=weight_vec[j]*((site_vec[j].t-meant)*(site_vec[j].t-meant)
                                    +(site_vec[j].u-meanu)*(site_vec[j].u-meanu)
                                    +(site_vec[j].v-meanv)*(site_vec[j].v-meanv)
                                    +(site_vec[j].w-meanw)*(site_vec[j].w-meanw)
                                    +(site_vec[j].x-meanx)*(site_vec[j].x-meanx)
                                    +(site_vec[j].y-meany)*(site_vec[j].y-meany)
                                    +(site_vec[j].z-meanz)*(site_vec[j].z-meanz));
		}
		fout<<s<<"\t"<<site_vec.size()<<"\t"<<rxy<<"\t"<<sqrt(rds_gyr/sumweight)<<std::endl;
	}
	fout.close();

	if(isotropic)
    {
        suffix=std::string("_sr.dat");
    }
    else
    {
        suffix=std::string("_anti_sr.dat");
    }
	ofname=std::string("7_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	fout.open(ofname.c_str());
	for(size_t j=0; j<sco; ++j)
    {
        if(n_array[j]>0)
        {
            fout<<j+1<<"\t"<<r_array[j]/n_array[j]<<std::endl;
        }
    }
	fout.close();
	delete[] r_array;
	delete[] n_array;
}

/********************* seven dimensional Bak-Sneppen model ***********************/
class site8d
{
public:
    int s;
    int t;
    int u;
    int v;
    int w;
    int x;
    int y;
    int z;

    site8d():s(0),t(0),u(0),v(0),w(0),x(0),y(0),z(0){}
    site8d(int _s, int _t, int _u, int _v, int _w, int _x, int _y, int _z) : s(_s), t(_t), u(_u), v(_v), w(_w), x(_x), y(_y), z(_z){}

    bool operator==(site8d const& other) const
    {
        return s==other.s && t==other.t && u==other.u && v==other.v && w==other.w && x==other.x && y==other.y && z==other.z;
    }

    friend std::size_t hash_value(site8d const& site)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, site.s);
        boost::hash_combine(seed, site.t);
        boost::hash_combine(seed, site.u);
        boost::hash_combine(seed, site.v);
        boost::hash_combine(seed, site.w);
        boost::hash_combine(seed, site.x);
        boost::hash_combine(seed, site.y);
        boost::hash_combine(seed, site.z);
        return seed;
    }
};

void bs8d(double f0, double alpha, int nei, int sco, int N, bool isotropic=true)
{
    MinHeap<double> h; // the heap to store sites
    std::vector<site8d> site_vec; // from handle to site, site_vec[handle]
    std::vector<int> weight_vec; // weight for each sites, i.e. covered times
    double *r_array=new double[sco]; // positive distance r relative to the origin
    int *n_array=new int[sco];
    boost::unordered_map<site8d, size_t> handle_map; // a map to from node to handle
    boost::unordered_map<site8d, size_t>::iterator handle_it;
	double r; // a temporary random number
	double sums,sumt,sumu,sumv,sumw,sumx,sumy,sumz,sumweight; // temporary variable for calculate center of mass
	double means,meant,meanu,meanv,meanw,meanx,meany,meanz; // center of mass of a single avalanche
	double rds_gyr; // radius of gyration
    int s,rxy,tph; // avalanche size and top handle
    std::string suffix; // suffix of the output file

	// initial r_array and n_array
	for(int i=0; i<sco; ++i)
	{
	    r_array[i]=0.0;
	    n_array[i]=0;
	}

	// set the output file name
	if(isotropic)
    {
        suffix=std::string("_snrR.dat");
    }
    else
    {
        suffix=std::string("_anti_snrR.dat");
    }
	std::string ofname=std::string("8_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	std::ofstream fout(ofname.c_str());

    for(int i=0;i<N;i++){

        // initial
        site_vec.clear();
		handle_map.clear();
		h.clear();
		s=0;
		tph=0;
		site8d min_site(0,0,0,0,0,0,0,0), nei_site;
		rxy=0;

        // set the active site
        site_vec.push_back(min_site);
        weight_vec.push_back(1);
		handle_map.insert(std::make_pair(min_site, h.push(0.0)));

        // start the avalanche
		while( h.topv()<f0 && s<sco ){

			s++;
			h.update(tph,dis()); // update min site

			// record the instances of having activity at a positive distance r relative to the origin at update step s
			r_array[s-1]+=abs(min_site.s)+abs(min_site.t)+abs(min_site.u)+abs(min_site.v)+abs(min_site.w)+abs(min_site.x)+abs(min_site.y)+abs(min_site.z);
			n_array[s-1]++;

            // ***********update neighbors*************//
			// update right neighbors along dimension 1
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site8d(min_site.s+j,min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }


			// update right neighbors along dimension 2
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site8d(min_site.s,min_site.t+j,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 3
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site8d(min_site.s,min_site.t,min_site.u+j,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 4
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v+j,min_site.w,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 5
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v,min_site.w+j,min_site.x,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 6
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v,min_site.w,min_site.x+j,min_site.y,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 7
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y+j,min_site.z);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

            // update right neighbors along dimension 8
			for(int j=1; j<=nei; ++j)
            {
                if(ui()<=alpha){
                    nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z+j);
                    handle_it=handle_map.find(nei_site);
                    r=dis();
                    if(handle_it!=handle_map.end()){
                        h.update(handle_it->second,r);
                        weight_vec[handle_it->second]++;
                    }
                    else{
                        if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                        {
                            rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                        }
                        site_vec.push_back(nei_site);
                        weight_vec.push_back(1);
                        handle_map.insert(std::make_pair(nei_site, h.push(r)));
                    }
                }
            }

			if(isotropic)
            {
                // update left neighbors along dimension 1
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site8d(min_site.s-j,min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 2
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site8d(min_site.s,min_site.t-j,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 3
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site8d(min_site.s,min_site.t,min_site.u-j,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 4
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v-j,min_site.w,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 5
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v,min_site.w-j,min_site.x,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 6
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v,min_site.w,min_site.x-j,min_site.y,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 7
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y-j,min_site.z);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }

                // update left neighbors along dimension 8
                for(int j=1; j<=nei; ++j)
                {
                    if(ui()<=alpha){
                        nei_site=site8d(min_site.s,min_site.t,min_site.u,min_site.v,min_site.w,min_site.x,min_site.y,min_site.z-j);
                        handle_it=handle_map.find(nei_site);
                        r=dis();
                        if(handle_it!=handle_map.end()){
                            h.update(handle_it->second,r);
                            weight_vec[handle_it->second]++;
                        }
                        else{
                            if(rxy<abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z))
                            {
                                rxy=abs(nei_site.s)+abs(nei_site.t)+abs(nei_site.u)+abs(nei_site.v)+abs(nei_site.w)+abs(nei_site.x)+abs(nei_site.y)+abs(nei_site.z);
                            }
                            site_vec.push_back(nei_site);
                            weight_vec.push_back(1);
                            handle_map.insert(std::make_pair(nei_site, h.push(r)));
                        }
                    }
                }
            }

			// reset
            tph=h.toph();
            min_site=site_vec[tph];
		}

		sums=0; sumt=0; sumu=0; sumv=0; sumw=0; sumx=0; sumy=0; sumz=0; sumweight=0;
		for(size_t j=0; j<site_vec.size(); ++j)
        {
            sums+=site_vec[j].s*weight_vec[j];
            sumt+=site_vec[j].t*weight_vec[j];
            sumu+=site_vec[j].u*weight_vec[j];
            sumv+=site_vec[j].v*weight_vec[j];
            sumw+=site_vec[j].w*weight_vec[j];
            sumx+=site_vec[j].x*weight_vec[j];
            sumy+=site_vec[j].y*weight_vec[j];
            sumz+=site_vec[j].z*weight_vec[j];
            sumweight+=weight_vec[j];
        }
        means=sums/sumweight;
        meant=sumt/sumweight;
        meanu=sumu/sumweight;
        meanv=sumv/sumweight;
        meanw=sumw/sumweight;
		meanx=sumx/sumweight;
		meany=sumy/sumweight;
		meanz=sumz/sumweight;
		rds_gyr=0;
		for(size_t j=0; j<site_vec.size(); ++j){
            rds_gyr+=weight_vec[j]*((site_vec[j].s-means)*(site_vec[j].s-means)
                                    +(site_vec[j].t-meant)*(site_vec[j].t-meant)
                                    +(site_vec[j].u-meanu)*(site_vec[j].u-meanu)
                                    +(site_vec[j].v-meanv)*(site_vec[j].v-meanv)
                                    +(site_vec[j].w-meanw)*(site_vec[j].w-meanw)
                                    +(site_vec[j].x-meanx)*(site_vec[j].x-meanx)
                                    +(site_vec[j].y-meany)*(site_vec[j].y-meany)
                                    +(site_vec[j].z-meanz)*(site_vec[j].z-meanz));
		}
		fout<<s<<"\t"<<site_vec.size()<<"\t"<<rxy<<"\t"<<sqrt(rds_gyr/sumweight)<<std::endl;
	}
	fout.close();

	if(isotropic)
    {
        suffix=std::string("_sr.dat");
    }
    else
    {
        suffix=std::string("_anti_sr.dat");
    }
	ofname=std::string("8_")+num2str(nei)+std::string("_")+num2str(alpha)+std::string("_")+num2str(f0)+suffix;
	fout.open(ofname.c_str());
	for(size_t j=0; j<sco; ++j)
    {
        if(n_array[j]>0)
        {
            fout<<j+1<<"\t"<<r_array[j]/n_array[j]<<std::endl;
        }
    }
	fout.close();
	delete[] r_array;
	delete[] n_array;
}


#endif // BS_H
