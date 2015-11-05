#include <globals.h>
#include <e_compon.h>

namespace {
using std::string;

class DEV_G_POLY_K : public COMPONENT { //
	public:
		string value_name() const { return "dummy"; }
		string port_name(uint_t) const { return "dummy"; }

		CARD* clone()const
		{ untested();
			const CARD* c = device_dispatcher["cpoly_g"];
			assert(c);
			CARD* c2 = c->clone();
			COMPONENT* d = prechecked_cast<COMPONENT*>(c2);
			assert(d);
			const COMMON_COMPONENT* b = bm_dispatcher["poly_k"];
			assert(b);
			COMMON_COMPONENT* bc = b->clone();
			d->attach_common(bc);
			// d->set_dev_type("g_poly_k");
			// assert(d->dev_type() == "g_poly_k");
			return d;
		}
}p1;

DISPATCHER<CARD>::INSTALL d1(&device_dispatcher, "g_poly_k|G_poly", &p1);
}
