#ifndef __MPIUTILITIES
#define __MPIUTILITIES

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <sstream>
#include <json_spirit.h>

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
	
	void read_json(std::string name, json_spirit::mValue& temp) {
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
		
		json_spirit::read(buffer, temp);
	}	
	
	void write_json(json_spirit::mValue const& temp, std::string name) {
		if(rank() == master) {
			std::ofstream file(name.c_str());
		    json_spirit::write(temp, file, json_spirit::pretty_print | json_spirit::single_line_arrays | json_spirit::remove_trailing_zeros);
			file.close();
		}
#ifdef HAVE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	};



	void read_json_all_processors(std::string name, json_spirit::mValue& temp) {
		int size;
		std::string buffer;
		
		std::ifstream file(name.c_str()); 
		if(!file) throw std::runtime_error(name + " not found.");
		
		file.seekg(0, std::ios::end);
		size = file.tellg();
		file.seekg(0, std::ios::beg);
		
		buffer.resize(size);
		file.read(&buffer[0], size);
		
		file.close();	
		
		json_spirit::read(buffer, temp);
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
