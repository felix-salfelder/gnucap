// do nothing but define bogus version

static const unsigned iv = 42;
static const char* in = "bogus///";
extern "C" {
	const unsigned interface_version(){return iv;}
	const char* interface_name(){return in;}
};

static const int aaaaa=5;
