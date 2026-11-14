// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef DownstreamPropagator_h
#define DownstreamPropagator_h
#include "classes/DelphesModule.h"
#include "modules/DownstreamTransport.h"
#include "modules/BirthKinematics.h"
class TIterator;
class TObjArray;

// Truth-level forward dipole transport and SciFi envelope acceptance.
// No UT requirement, detector resolution, LLP decay generation or event selection.
class DownstreamPropagator: public DelphesModule {
public:
  DownstreamPropagator();
  ~DownstreamPropagator();
  void Init();
  void Process();
  void Finish();
private:
  downstream::Config fConfig; //!
  downstream::BirthCuts fBirthCuts; //!
  TIterator *fIterator; //!
  const TObjArray *fInputArray; //!
  TObjArray *fOutputArray; //!
  ClassDef(DownstreamPropagator, 1)
};
#endif
