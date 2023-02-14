// Karl Yerkes
// 2023-02-07
// MAT 240B ~ Audio Programming
//

#include <juce_audio_processors/juce_audio_processors.h>

#include <vector>

template <typename T>
T mtof(T m) {
  return T(440) * pow(T(2), (m - T(69)) / T(12));
}

template <typename T>
T dbtoa(T db) {
  return pow(T(10), db / T(20));
}

struct CircularArray {
  std::vector<float> _data;
  // index [0.0, size]
  float operator[](float index) {
    // safety
    index = wrap(index, _data.size(), 0);
    // need i, j, and t ~ given index
    int i = (int)index;
    int j = (i + 1 > _data.size() - 1) ? 0 : i + 1;
    float t = index - i;
    return _data[i] * (1 - t) + t * _data[i + 1];
  }
};

class DelayLine : std::vector<float> {
  //
  int index = 0;

 public:
  float read(float seconds_ago, float samplerate) {
    //
    jassert(seconds_ago < size() / samplerate);

    float i = index - seconds_ago * samplerate;
    if (i < 0) {
      i += size();
    }
    return at((int)i);  // no linear interpolation
  }

  void write(float value) {
    jassert(size() > 0);
    at(index) = value;  // overwrite the oldest value

    // handle the wrapping for circular buffer
    index++;
    if (index >= size()) index = 0;
  }

  void allocate(float seconds, float samplerate) {
    // floor(seconds * samplerate) + 1 samples
    resize((int)floor(seconds * samplerate) + 1);
  }
};

class Delay : public juce::AudioProcessor {
  juce::AudioParameterFloat* gain;
  juce::AudioParameterFloat* delay_time;
  juce::AudioParameterChoice* choice;

  /// add parameters here ///////////////////////////////////////////////////

  DelayLine delay_line;

 public:
  Delay()
      : AudioProcessor(
            BusesProperties()
                .withInput("Input", juce::AudioChannelSet::stereo())
                .withOutput("Output", juce::AudioChannelSet::stereo())) {
    addParameter(gain = new juce::AudioParameterFloat(
                     {"gain", 1}, "Gain",
                     juce::NormalisableRange<float>(-65, -1, 0.01f), -65));
    addParameter(delay_time = new juce::AudioParameterFloat(
                     {"delay", 1}, "Delay",
                     juce::NormalisableRange<float>(0, 4, 0.01f), 0.7f));
    addParameter(choice = new juce::AudioParameterChoice(
                     "choice", "Choice", juce::StringArray({"Foo", "Bar"}), 0));
    /// add parameters here /////////////////////////////////////////////

    // XXX getSampleRate() is not valid here

    // delay_line.allocate(delay_time->getNormalisableRange().end,
    //                     (float)getSampleRate());

    // std::cout << delay_time->getNormalisableRange().end << std::endl;
  }

  /// handling the actual audio! ////////////////////////////////////////////
  void processBlock(juce::AudioBuffer<float>& buffer,
                    juce::MidiBuffer&) override {
    buffer.clear(0, 0, buffer.getNumSamples());
    auto left = buffer.getWritePointer(0, 0);
    auto right = buffer.getWritePointer(1, 0);

    //   printf("got here\n"); // i/o ~ might take a while

    for (int i = 0; i < buffer.getNumSamples(); ++i) {
      delay_line.write((left[i] + right[i]) / 2);
      left[i] = delay_line.read(delay_time->get(), (float)getSampleRate());
      right[i] = left[i];
    }
  }

  /// handle doubles ? //////////////////////////////////////////////////////
  // void processBlock(AudioBuffer<double>& buffer, MidiBuffer&) override {
  //   buffer.applyGain(dbtoa((float)*gain));
  // }

  /// start and shutdown callbacks///////////////////////////////////////////
  void prepareToPlay(double samplerate, int) override {
    // XXX when does this get called? seems to not get called in stand-alone
    delay_line.allocate(delay_time->getNormalisableRange().end,
                        (float)samplerate);
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
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Delay)
};

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() { return new Delay(); }
