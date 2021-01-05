
#include <multem.h>

int main() {
  
  mt::SystemConfiguration a;
  mt::Input<float> b;
  b.set_system_conf(a);
  mt::test(b);

  return 0;
}
