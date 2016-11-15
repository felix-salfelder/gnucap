#include <globals.h>
#include <e_model.h>

class CARD;

namespace {
class TEST : public MODEL_CARD{
public:
  TEST() : MODEL_CARD(NULL){};
  CARD*	 clone()const{ return new TEST(*this);};
  std::string value_name()const{ return ""; }
} a;

DISPATCHER<MODEL_CARD>::INSTALL d1(&model_dispatcher, "Pmos8", &a);
}

#include <stdio.h>
