#ifndef __MPIUTILITIES
#define __MPIUTILITIES

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include <valarray>
#include "IO.h"
using json=nlohmann::json;

namespace mpi {		
	//-----------------------------------------------------------------------------------------------------------------------------------------------------
	int const master = 0;
	
	int rank() {
		int temp = 0;
		
#ifdef HAVE_MPI		
		MPI_Comm_rank(MPI_COMM_WORLD, &temp);
#endif
		
		return temp;
	};
	
	int number_of_workers() {
		int temp = 1;
		
#ifdef HAVE_MPI
		MPI_Comm_size(MPI_COMM_WORLD, &temp);
#endif
		
		return temp;
	};
	
	void read_json(std::string name, json& jObject) {
		int size;
		std::string buffer;
		
		if(rank() == master) {
			std::ifstream file(name.c_str()); 
			if(!file) throw std::runtime_error(name + " not found.");
			
			file.seekg(0, std::ios::end);
			size = file.tellg();
			file.seekg(0, std::ios::beg);
			
			buffer.resize(size);
			file.read(&buffer[0], size);
			
			file.close();
		}			
		
#ifdef HAVE_MPI
		MPI_Bcast(&size, 1, MPI_INT, master, MPI_COMM_WORLD);
		
		if(rank() != master) buffer.resize(size);
		
		MPI_Bcast(&buffer[0], size, MPI_CHAR, master, MPI_COMM_WORLD); 
#endif
		
		jObject = json::parse(&buffer[0],&buffer[0] + size);
	}	
	
	void write_json(std::string name,json const& jObject) {
		if(rank() == master) {
			IO::writeJsonToFile(name,jObject,true);
		}
#ifdef HAVE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	};
	void reduce(std::valarray<double> &toTransmit,std::valarray<double> &toReceive){
#ifdef HAVE_MPI
			MPI_Reduce(&toTransmit[0], mpi::rank() == mpi::master ? &toReceive[0] : 0, toTransmit.size(), MPI_DOUBLE, MPI_SUM, mpi::master, MPI_COMM_WORLD);	
#else
			toReceive = toTransmit;
#endif
	}
	void reduce(uint64_t &toTransmit,uint64_t &toReceive){
#ifdef HAVE_MPI
			MPI_Reduce(&toTransmit, mpi::rank() == mpi::master ? &toReceive : 0, 1, MPI_DOUBLE, MPI_SUM, mpi::master, MPI_COMM_WORLD);	
#else
			toReceive = toTransmit;
#endif
	}
	void allReduce(std::valarray<double> &toTransmit,std::valarray<double> &toReceive){
#ifdef HAVE_MPI
			MPI_Allreduce(&toTransmit[0], &toReceive[0], toTransmit.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	
#else
			toReceive = toTransmit;
#endif
	}
	void allReduce(uint64_t &toTransmit,uint64_t &toReceive){
#ifdef HAVE_MPI
			MPI_Allreduce(&toTransmit, &toReceive, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	
#else
			toReceive = toTransmit;
#endif
	}

	void getMax(uint64_t &toTransmit,uint64_t &toReceive){
#ifdef HAVE_MPI
			MPI_Allreduce(&toTransmit, &toReceive, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);	
#else
			toReceive = toTransmit;
#endif
	}


	void read_json_all_processors(std::string name, json& jObject) {
		IO::readJsonFile(name,jObject);
	}	
	
// This is a temporary solution !!!!!!
	struct ofstream : public std::ostringstream {
		explicit ofstream(std::string filename, std::ios_base::openmode mode = ios_base::out) : std::ostringstream(mode), filename_(filename), mode_(mode) {
		};
		~ofstream() {
			if(rank() == master) {
			    std::ofstream file(filename_.c_str(), mode_);
				if(!file)
					throw std::runtime_error("mpi::ofstream: Can not open " + filename_ + " file !");
				file << str();
				file.close();
			}
#ifdef HAVE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
		}
	private:
		std::string filename_;
		ios_base::openmode mode_;
	};
	
	enum State { every, one, none };
	
	struct Cout {
		Cout() : state_(one) {};
		void operator=(State state) {state_ = state;};
		std::ostream& operator<<(std::ostream& (*pf)(std::ostream&)) { 
			if(state_ == every) return pf(std::cout << rank() << " : "); 
			if(state_ == one) return mpi::rank() == mpi::master ? pf(std::cout) : null_;
			return null_;
		};
	private:
		static std::ostream null_;
		State state_;
		
		template<typename T> 
		friend std::ostream& operator<<(Cout const&, T);
	};
	
	std::ostream Cout::null_(0);
	
	template<typename T>
	std::ostream& operator<<(Cout const& c, T t) {
		if(c.state_ == every) return std::cout << rank() << " : " << t;
		if(c.state_ == one) return mpi::rank() == mpi::master ? std::cout << t : Cout::null_;
		return Cout::null_;
	};
	
	extern Cout cout;
	Cout cout;
};


#endif
