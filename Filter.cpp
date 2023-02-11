// Karl Yerkes
// 2023-02-07
// MAT 240B ~ Audio Programming
//

#include <juce_audio_processors/juce_audio_processors.h>

template <typename T>
T mtof(T m) {
  return T(440) * pow(T(2), (m - T(69)) / T(12));
}

template <typename T>
T dbtoa(T db) {
  return pow(T(10), db / T(20));
}

class BiquadFilter {
  // Audio EQ Cookbook
  // http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt

  // x[n-1], x[n-2], y[n-1], y[n-2]
  float x1 = 0, x2 = 0, y1 = 0, y2 = 0;

  // filter coefficients
  float b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0;

 public:
  float operator()(float x0) {
    // Direct Form 1, normalized...
    float y0 = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
    y2 = y1;
    y1 = y0;
    x2 = x1;
    x1 = x0;
    return y0;
  }

  void normalize(float a0) {
    b0 /= a0;
    b1 /= a0;
    b2 /= a0;
    a1 /= a0;
    a2 /= a0;
  }

  void lpf(float f0, float Q, float samplerate) {
    float w0 = 2 * float(M_PI) * f0 / samplerate;
    float alpha = sin(w0) / (2 * Q);
    b0 = (1 - cos(w0)) / 2;
    b1 = 1 - cos(w0);
    b2 = (1 - cos(w0)) / 2;
    float a0 = 1 + alpha;
    a1 = -2 * cos(w0);
    a2 = 1 - alpha;

    normalize(a0);
  }
};

class History {
  float _data = 0;

 public:
  float operator()(float in) {
    float value = _data;
    _data = in;
    return value;
  }
  float operator()() { return _data; }
};

// adapted from this JOS3 paper:
// https://ccrma.stanford.edu/~jos/svf/svf.pdf
//
class StateVariableFilter {
  // TODO:
  // * parameterize q and wct
  // * compress expressions
  // * rename variables
  History z1;
  History z2;

  float _high = 0;
  float _mid = 0;
  float _low = 0;

 public:
 void zero() { _data = 0; }
  void step(float in, float wct, float q) {
    _mid = z2();
    float mul_3 = z2() * -q;
    float mul_4 = z2() * wct;
    float add_5 = mul_4 + z1();
    _low = add_5;
    float mul_6 = add_5 * -1;
    float add_7 = mul_3 + mul_6;
    float add_8 = in + add_7;
    _high = add_8;
    float mul_9 = add_8 * wct;
    float add_10 = mul_9 + z2();
    // https://stackoverflow.com/questions/2487653/avoiding-denormal-values-in-c
    // float history_1_next_11 = fixdenorm(add_5);
    // float history_2_next_12 = fixdenorm(add_10);
    // float z1 = history_1_next_11;
    // float z2 = history_2_next_12;
    z1(add_5);
    z2(add_10);
    // for a clearer implementation, go here:
    // https://github.com/JordanTHarris/VAStateVariableFilter/blob/master/Source/Effects/VAStateVariableFilter.cpp
    // of maybe look at this:
    // https://github.com/JamesWenlock/StateVariableFilter/blob/master/Source/PluginProcessor.cpp
  }

  float high() { return _high; }
  float mid() { return _mid; }
  float low() { return _low; }

  // addapted from this Max/Gen implementation of the JOS3 paper
  /*
<pre><code>
----------begin_max5_patcher----------
646.3ocyW1saaCBFF9X6qBKNbKMxefM1d2JSUStInVpRvV1jtTU068Y9wSoU
taPMTkSZEDBu7v2O7lWRSP20clMhx9Q1OyRRdIMIQOkZhD63Dzw1y6NzNpWF
Rv9c2cOh1X9HI6rTOMWjQlmTb5X2I4AlT+M.6rlojO2yLxgPY2Z+H9d8dLsu
2.UWrKbw7ljamruUt6At39eMv1IM6SMTuMeSFPHp+g0Cv3s4Y2p9Jullp9yl
0Q22BBZzEQC+wnUlCJZJMDha1VdshVounUPKUzPxqhGZSm8L7x3ku.CEKx.7
OBOkkpCukAR40alGw2vCEnWfV7BOD2COXeCOzFcYiJSKdgmraffDg.uoyDgJ
v5VD3pXP22CBZ4d21qH+h.WSQLR9BBZMd20yjENSFDCxdfOJ6FdND7U8I6Hh
wFLKtZSJq8+o3pKHKNMDCXji5c6DpNXQpfq7HW4m8oLCYwIxM49MLOC3q4Wa
CkY3hh42I3vg.NruvUUnoBr0bQANkIDvYSHdaQzFerN3+.OH5KKzAt38+RL8
Fpl+sbM1cZX27UsslHK+uW46YiRtnUx6DWrF5aVS2vd1f9ru3UoqBScPX0gC
9fPlq5T4fNP4xDRVkx0Nnr8zEcBqM2juGPb7uZKhQxiSJSVV47UoL3b9zZSb
AWJMUNjWcBDPcLCZsB03BQzXD1vtpLD5TUWXNJ0HXWqQffWi3fxMgHuk33aH
qUHmpDwQocC3ZS.3qfQXYFW2yHy8v9u4KvBVeZ66ehMLZWrViIyaO1oOW0az
C4ByvR8vA1S740q89iZGl7dImLdcZPetPmoF2ZnicSDJNwsUiSzMIo1Xnn8H
aru0.BZ+X+16YBT5qo+A9hsJf.
-----------end_max5_patcher-----------
</code></pre>
  */
};

// using namespace juce;

class KarplusStrong : public juce::AudioProcessor {
  juce::AudioParameterFloat* gain;
  juce::AudioParameterFloat* note;
  juce::AudioParameterFloat* q;
  /// add parameters here ///////////////////////////////////////////////////

  // BiquadFilter filter;
  StateVariableFilter filter;

 public:
  KarplusStrong()
      : AudioProcessor(
            BusesProperties()
                .withInput("Input", juce::AudioChannelSet::stereo())
                .withOutput("Output", juce::AudioChannelSet::stereo())) {
    addParameter(gain = new juce::AudioParameterFloat(
                     {"gain", 1}, "Gain",
                     juce::NormalisableRange<float>(-65, -1, 0.01f), -65));
    addParameter(note = new juce::AudioParameterFloat(
                     {"note", 1}, "Note",
                     juce::NormalisableRange<float>(-2, 129, 0.01f), 40));
    addParameter(
        q = new juce::AudioParameterFloat(
            {"q", 1}, "Q", juce::NormalisableRange<float>(0, 4, 0.01f), 0.7f));
    /// add parameters here /////////////////////////////////////////////

    // XXX getSampleRate() is not valid here
  }

  /// handling the actual audio! ////////////////////////////////////////////
  void processBlock(juce::AudioBuffer<float>& buffer,
                    juce::MidiBuffer&) override {
    buffer.clear(0, 0, buffer.getNumSamples());
    auto left = buffer.getWritePointer(0, 0);
    auto right = buffer.getWritePointer(1, 0);

    for (int i = 0; i < buffer.getNumSamples(); ++i) {
      // filter.lpf(mtof(note->get()), q->get(), (float)getSampleRate());
      // left[i] = filter((left[i] + right[i]) / 2) * dbtoa(gain->get());
      filter.step((left[i] + right[i]) / 2, note->get() / 127, q->get());
      if (std::isnan(f)) {
        filter.reset();
        f = filter.low();
      }
      left[i] = filter.low();
      right[i] = left[i];
    }
  }

  /// handle doubles ? //////////////////////////////////////////////////////
  // void processBlock(AudioBuffer<double>& buffer, MidiBuffer&) override {
  //   buffer.applyGain(dbtoa((float)*gain));
  // }

  /// start and shutdown callbacks///////////////////////////////////////////
  void prepareToPlay(double, int) override {
    // XXX when does this get called? seems to not get called in stand-alone
  }
  void releaseResources() override {}

  /// maintaining persistant state on suspend ///////////////////////////////
  void getStateInformation(juce::MemoryBlock& destData) override {
    juce::MemoryOutputStream(destData, true).writeFloat(*gain);
  }

  void setStateInformation(const void* data, int sizeInBytes) override {
    gain->setValueNotifyingHost(
        juce::MemoryInputStream(data, static_cast<size_t>(sizeInBytes), false)
            .readFloat());
  }

  /// general configuration /////////////////////////////////////////////////
  const juce::String getName() const override { return "Quasi Band Limited"; }
  double getTailLengthSeconds() const override { return 0; }
  bool acceptsMidi() const override { return true; }
  bool producesMidi() const override { return false; }

  /// for handling presets //////////////////////////////////////////////////
  int getNumPrograms() override { return 1; }
  int getCurrentProgram() override { return 0; }
  void setCurrentProgram(int) override {}
  const juce::String getProgramName(int) override { return "None"; }
  void changeProgramName(int, const juce::String&) override {}

  /// ?????? ////////////////////////////////////////////////////////////////
  bool isBusesLayoutSupported(const BusesLayout& layouts) const override {
    const auto& mainInLayout = layouts.getChannelSet(true, 0);
    const auto& mainOutLayout = layouts.getChannelSet(false, 0);

    return (mainInLayout == mainOutLayout && (!mainInLayout.isDisabled()));
  }

  /// automagic user interface //////////////////////////////////////////////
  juce::AudioProcessorEditor* createEditor() override {
    return new juce::GenericAudioProcessorEditor(*this);
  }
  bool hasEditor() const override { return true; }

 private:
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(KarplusStrong)
};

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
  return new KarplusStrong();
}
