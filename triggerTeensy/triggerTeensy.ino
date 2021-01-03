/* #################################################
  ############## PIN CONFIGURATION ###################
  #################################################### */
// TTL outputs
// pins to create triggers
#define PIN_TRIALTRIG 11 // trial-start trigger that can be switched by serial command 'MAKE_TRIALTRIGGER'
#define PIN_STIMTRIG 12 // stimulus trigger that can be switched by serial command 'MAKE_STIMTRIGGER'
#define PIN_CAMTRIG 13 // camera trigger that can be switched by serial command 'MAKE_CAMTRIGGER'

// output for FISBA module 1
#define PIN_GND_1 22 // make this the ground pin for laser 1
#define PIN_RED_1 21 // trigger that is used to enable red (625nm) output from laser 1
#define PIN_BLUE_1 20 // trigger that is used to enable violet (405nm) output from laser 1 - CARFEUL this is referred to as 'blue' in the FISBA nomenclature but really violet
#define PIN_CYAN_1 19 // trigger that is used to enable cyan (488nm) output from laser 1 - CARFEUL this is referred to as 'green' in the FISBA nomenclature but its really more blue-ish and should be used for most optogenetics

// output for FISBA module 2
#define PIN_GND_2 17 // make this the ground pin for laser 2
#define PIN_RED_2 16 // trigger that is used to enable red (625nm) output from laser 2
#define PIN_BLUE_2 15 // trigger that is used to enable violet (405nm) output from laser 2 - CARFEUL this is referred to as 'blue' in the FISBA nomenclature but really violet
#define PIN_CYAN_2 14 // trigger that is used to enable cyan (488nm) output from laser 2 - CARFEUL this is referred to as 'green' in the FISBA nomenclature but its really more blue-ish and should be used for most optogenetics


// TTL inputs
// pins to read potential enable triggers from
#define PIN_ENABLE_1 8 // enable signal for laser 1
#define PIN_ENABLE_2 9 // enable signal for laser 2

/* #################################################
  ############## SERIAL COMMANDS ###################
  #################################################### */
// Byte codes for serial communication
// inputs
#define MAKE_STIMTRIGGER 101 // identifier to produce a stimulus trigger
#define MAKE_TRIALTRIGGER 102 // identifier to produce a trial onset trigger
#define MAKE_CAMTRIGGER 103 // identifier to produce a trigger for behavioral video cameras
#define CHANGE_ENABLETRIGGERS 150 // identifier to adjust trigger output lines (expects two subsequent bytes: a number between 1-3 to enable either red (1), cyan(2) or violet(3) output for each module

#define GOT_BYTE 10 // positive handshake for module_info
#define GOT_TRIGGER 11 // positive handshake for trigger command
#define DID_ABORT 15 // negative handshake for incoming bytes

// Other variables
bool midRead = false;
int FSMheader = 0;
unsigned long clocker = millis();
unsigned long camClocker = millis();
unsigned long trialClocker = millis();
unsigned long stimClocker = millis();
int stimDur = 10; // duration of stimulus trigger in ms
int trialDur = 50; // duration of trial trigger in ms
int camDur = 200; // duration of camera trigger in ms
bool stimTrigger = false;
bool trialTrigger = false;
bool camTrigger = false;
byte cByte = 0; // temporary variable for serial communication
volatile int enable_1 = PIN_RED_1; //current enable line for laser module 1
volatile int enable_2 = PIN_RED_2; //current enable line for laser module 2

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setup() {

  Serial.begin(11520);

  // Set pin modes for digital output lines
  pinMode(PIN_STIMTRIG, OUTPUT);
  pinMode(PIN_TRIALTRIG, OUTPUT);
  pinMode(PIN_CAMTRIG, OUTPUT);

  pinMode(PIN_GND_1, OUTPUT);
  pinMode(PIN_GND_2, OUTPUT);
  digitalWriteFast(PIN_GND_1, LOW);
  digitalWriteFast(PIN_GND_2, LOW);

  pinMode(PIN_RED_1, OUTPUT);
  pinMode(PIN_BLUE_1, OUTPUT);
  pinMode(PIN_CYAN_1, OUTPUT);
  pinMode(PIN_RED_2, OUTPUT);
  pinMode(PIN_BLUE_2, OUTPUT);
  pinMode(PIN_CYAN_2, OUTPUT); 

  // Set pin modes for digital input lines
  pinMode(PIN_ENABLE_1, INPUT);
  pinMode(PIN_ENABLE_2, INPUT);

  attachInterrupt(digitalPinToInterrupt(PIN_ENABLE_1), EnableChange_1, CHANGE); // interupt to create enable triggers
  attachInterrupt(digitalPinToInterrupt(PIN_ENABLE_2), EnableChange_2, CHANGE); // interupt to create enable triggers

}

void loop() {
  if (Serial.available() > 0) {
    if (!midRead) {
      FSMheader = Serial.read();
      midRead = true; // flag for current reading of serial information
      clocker = millis(); // counter to make sure that all serial information arrives within a reasonable time frame (currently 100ms)
    }

    if (FSMheader == MAKE_STIMTRIGGER) { // create stimulus trigger
      stimTrigger = true;
      stimClocker = millis();
      digitalWriteFast(PIN_STIMTRIG, HIGH); // set stimulus trigger to high

      Serial.write(GOT_TRIGGER);
      midRead = false;
    }

    if (FSMheader == MAKE_TRIALTRIGGER) { // create trial-onset trigger
      trialTrigger = true;
      trialClocker = millis();
      digitalWriteFast(PIN_TRIALTRIG, HIGH); // set trial trigger to high

      Serial.write(GOT_TRIGGER);
      midRead = false;
    }
    
    if (FSMheader == MAKE_CAMTRIGGER) { // create trial-onset trigger
      camTrigger = true;
      camClocker = millis();
      digitalWriteFast(PIN_CAMTRIG, HIGH); // set camera trigger to high

      Serial.write(GOT_TRIGGER);
      midRead = false;
    }

    if (FSMheader == CHANGE_ENABLETRIGGERS) { // check which enable lines should be used for laser modules 1 and 2

      if (Serial.available() > 1){
        // check laser l
        cByte = Serial.read();
        if (cByte == 1) {enable_1 = PIN_RED_1;} // enable red line
        else if (cByte == 2) {enable_1 = PIN_CYAN_1;} // enable 488nm line
        else if (cByte == 3) {enable_1 = PIN_BLUE_1;} // enable 405nm line
        else {Serial.write(DID_ABORT);}

        // check laser 2
        cByte = Serial.read();
        if (cByte == 1) {enable_2 = PIN_RED_2;} // enable red line
        else if (cByte == 2) {enable_2 = PIN_CYAN_2;} // enable 488nm line
        else if (cByte == 3) {enable_2 = PIN_BLUE_2;} // enable 405nm line
        else {Serial.write(DID_ABORT);}

        // done
        Serial.write(GOT_TRIGGER);
        midRead = false;
      }
    }
  }

  if (midRead && ((millis() - clocker) >= 1000)) {
    midRead = false; Serial.write(FSMheader);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // check stim trigger
  if (stimTrigger) {
    if ((millis() - stimClocker) > stimDur) {  // done with stim trigger
      digitalWriteFast(PIN_STIMTRIG, LOW); // set stimulus trigger to low
      stimTrigger = false;
    }
  }

  // check trial trigger
  if (trialTrigger) {
    if ((millis() - trialClocker) > trialDur) {  // done with trial trigger
      digitalWriteFast(PIN_TRIALTRIG, LOW); // set trial trigger to low
      trialTrigger = false;
    }
  }
  
  // check camera trigger
  if (camTrigger) {
    if ((millis() - camClocker) > camDur) {  // done with camera trigger
      digitalWriteFast(PIN_CAMTRIG, LOW); // set trial trigger to low
      camTrigger = false;
    }
  }
}

// interrupts
void EnableChange_1() {
  if (digitalReadFast(PIN_ENABLE_1)) {
    digitalWriteFast(enable_1, HIGH);
  }
  else {
      digitalWriteFast(enable_1, LOW);
  }
}

void EnableChange_2() {
  if (digitalReadFast(PIN_ENABLE_2)) {
    digitalWriteFast(enable_2, HIGH);
  }
  else {
      digitalWriteFast(enable_2, LOW);
  }
}  
