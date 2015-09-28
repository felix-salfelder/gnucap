/*--------------------------------------------------------------------------
 * expect makes the following simulation commands expect output
 * this is for testing purposes. see ../tests
 * (does not really make sense currently)
 */
class CMD_EXPECT : public CMD {
public:
  virtual ~CMD_EXPECT(){

  }
  void do_it(CS& cmd, CARD_LIST* )
  {
    trace0("CMD_EXPECT::do_it");
    unsigned here = cmd.cursor();
    try {
      std::string file_name;
      cmd >> file_name;
      CS* file = new CS(CS::_INC_FILE, file_name);
      trace1( (" CMD_EXPECT::do_it > " +file_name).c_str() , (long int)(OPT::language) );

      _sim->expect(file);


    }catch (Exception_File_Open& e) {
      cmd.warn(bDANGER, here, e.message() + '\n');
    }catch (Exception_End_Of_Input& e) {
      // done
    }
  }
} p0x;
DISPATCHER<CMD>::INSTALL d0x(&command_dispatcher, "expect", &p0x);
