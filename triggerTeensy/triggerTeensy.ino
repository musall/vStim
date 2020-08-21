/* #################################################
############## PIN CONFIGURATION ###################
#################################################### */
// TTL Outputs
#define PIN_TRIALTRIG 11 // trial-start trigger that can be switched by serial command 'MAKE_TRIALTRIGGER'
#define PIN_STIMTRIG 12 // stimulus trigger that can be switched by serial command 'MAKE_STIMTRIGGER'

// Byte codes for serial communication
// inputs
#define MODULE_INFO 50 // byte to return handshake
#define MAKE_STIMTRIGGER 101 // identifier to produce a stimulus trigger
#define MAKE_TRIALTRIGGER 102 // identifier to produce a trial onset trigger

#define GOT_BYTE 10 // positive handshake for module_info
#define GOT_TRIGGER 11 // positive handshake for trigger command
#define DID_ABORT 15 // negative handshake for incoming bytes

// Other variables
bool midRead = false;
int FSMheader = 0;
unsigned long clocker = millis();
unsigned long trialClocker = millis();
unsigned long stimClocker = millis();
int stimDur = 10; // duration of stimulus trigger in ms
int trialDur = 50; // duration of trial trigger in ms
bool stimTrigger = false;
bool trialTrigger = false;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setup() {
  
  Serial.begin(11520);
  
  // Set pin modes for digital output lines
  pinMode(PIN_STIMTRIG, OUTPUT);
  pinMode(PIN_TRIALTRIG, OUTPUT);

}

void loop() {
  if (Serial.available() > 0) {
    if (!midRead) {
      FSMheader = Serial.read();
      midRead = true; // flag for current reading of serial information
      clocker = millis(); // counter to make sure that all serial information arrives within a reasonable time frame (currently 100ms)
    }
   
    if (FSMheader == MODULE_INFO) { // return module information to bpod
      Serial.write(GOT_BYTE);
      midRead = false;
    }
    
    else if (FSMheader == MAKE_STIMTRIGGER) { // create stimulus trigger
      stimTrigger = true;
      stimClocker = millis();
      digitalWriteFast(PIN_STIMTRIG, HIGH); // set stimulus trigger to high

      Serial.write(GOT_TRIGGER);
      midRead = false;
    }
    
    else if (FSMheader == MAKE_TRIALTRIGGER) { // create trial-onset trigger
      trialTrigger = true;
      trialClocker = millis();
      digitalWriteFast(PIN_TRIALTRIG, HIGH); // set trial trigger to high

      Serial.write(GOT_TRIGGER);
      midRead = false;
    }
  }

  if (midRead && ((millis() - clocker) >= 100)) {
    midRead = 0; Serial.write(DID_ABORT);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // make stim trigger
  if (stimTrigger) {
    if ((millis() - stimClocker) > stimDur) {  // done with stim trigger
      digitalWriteFast(PIN_STIMTRIG, LOW); // set stimulus trigger to low
      stimTrigger = false;
    }
  }

  // make trial trigger
  if (trialTrigger) {
    if ((millis() - trialClocker) > trialDur) {  // done with trial trigger
      digitalWriteFast(PIN_TRIALTRIG, LOW); // set trial trigger to low
      trialTrigger = false;
    }
  }
}
