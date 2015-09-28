#include <globals.h>
#include <e_compon.h>

// this is a test only. bm_wrapper is now part of gnucap-bm

namespace {
using std::string;

class DEV_VS_SIN : public COMPONENT { //
	public:
		string value_name() const { return "dummy"; }
		string port_name(uint_t) const { return "dummy"; }

		CARD* clone()const
		{ untested();
			const CARD* c = device_dispatcher["V"];
			assert(c);
			CARD* c2 = c->clone();
			COMPONENT* d = prechecked_cast<COMPONENT*>(c2);
			assert(d);
			const COMMON_COMPONENT* b = bm_dispatcher["sin"];
			assert(b);
			COMMON_COMPONENT* bc = b->clone();
			d->attach_common(bc);
			d->set_dev_type("vsource_sin");
			assert(d->dev_type() == "vsource_sin");
			return d;
		}
}p1;

DISPATCHER<CARD>::INSTALL d1(&device_dispatcher, "vsource_sin", &p1);
}
