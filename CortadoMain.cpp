#include "./OpenGLTerminalManager.h"
#include "Cortado.h"

int main() {
  OpenGLTerminalManager terminalManager("Cortado");
  Cortado cortado(&terminalManager);
  cortado.run();
  return 0;
}
