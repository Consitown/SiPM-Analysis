{
	std::cout << "\n --- Loading wavecatcher analysis library --- \n\n";
	// load library only if system is not windows
	if (std::string(gSystem->GetBuildArch()).find("win") == std::string::npos) gSystem->Load("ReadRunLib.sl");
	std::cout << "\nIf you get an error about ReadRunLib.sl please compile the library again or adjust the path to ReadRunLib.sl in rootlogon.c\n\n";
}