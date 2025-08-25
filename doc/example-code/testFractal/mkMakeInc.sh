#!/bin/sh

print_usage_and_exit () {

    echo "$0 generates a make.inc file, which is used by the Makefile"
    echo "to determine the locations of the KMCThinFilm library and Boost,"
    echo "and the C++ compiler to use."
    echo
    echo "Available options are as follows:"
    echo "    --help                Prints this message and exits"
    echo "    --interactive or -i   Generates the make.inc file from questions"
    echo "                          asked at the command prompt [default]"
    echo "    --batch or -b         Generates a make.inc with dummy values for"
    echo "                          the locations of KMCThinFilm library, etc."
    echo "                          One may then edit the make.inc file afterwards."

    exit 1
}

choose_compiler_options () {

    CXX="$1"
    
    CXXFLAGS_GNU="-Wall -O3"
    CXXFLAGS_INTEL="-Wall -O3" # These happen to be the same options
			       # as for the Gnu compiler for now.

    if [ "x$CXX" = "xg++" ]
    then
	CXXFLAGS="$CXXFLAGS_GNU"
    elif [ "x$CXX" = "xicpc" ]
    then
	CXXFLAGS="$CXXFLAGS_INTEL"
    else
	CXXFLAGS=""
    fi

    echo "$CXXFLAGS"
}

# Dummy default values
KMC_INST="/path/to/desired/KMCThinFilm/install/location/inst"
BOOST_ROOT="/path/to/Boost"
CXX="g++"
CXXFLAGS=`choose_compiler_options "$CXX"`

INTERACTIVE=1

while [ "$#" -gt 0 ]
do
    case "$1" in
	--help)
	    print_usage_and_exit
	    ;;
	--interactive | -i)
	    INTERACTIVE=1
	    ;;
	--batch | -b)
	    INTERACTIVE=0
	    ;;
	*)
	    echo "Unrecognized argument: $1"
	    print_usage_and_exit
	    ;;
    esac
    shift
done

if [ $INTERACTIVE -eq 1 ]
then
    echo "What is the root of your installation of the KMCThinFilm library?"
    read -r KMC_INST

    echo "What is the root of your installation of Boost?"
    read -r BOOST_ROOT
    
    echo "What is your C++ compiler?"
    read -r CXX

    CXXFLAGS=`choose_compiler_options "$CXX"`

    if [ -z "$CXXFLAGS" ]
    then
	CXXFLAGS=unknown
    fi
    
    BAD_ANSWER=1 
    while [ $BAD_ANSWER -eq 1 ]
    do	
	echo "Default compiler option(s): $CXXFLAGS. Is this okay? Answer Y for yes and N for no."
	read -r ANSWER

	case "$ANSWER" in
	    y*|Y*)
		BAD_ANSWER=0
		;;
	    n*|N*)
		BAD_ANSWER=0
		echo "Which compiler options do you wish to use?"
		read -r CXXFLAGS
		;;
	    *)
		echo "Bad answer: $ANSWER; Answer Y for yes and N for no."
		;;
	esac
    done

fi

cat > make.inc <<EOF
# These paths should be set to the actual paths needed for compilation
KMC_INST = $KMC_INST
BOOST_ROOT = $BOOST_ROOT

# Set the compiler to one of those available at one's workstation or cluster.
CXX = $CXX

# Options for the compiler
CXXFLAGS = $CXXFLAGS
EOF
