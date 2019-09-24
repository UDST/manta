enum Phase {Initialization, Routing, Simulation, Export};

class DataExporter
{

public:
    // TODO: Start in Initialization mode
    void DataExporter(void);

    // TODO: Assert current_phase_ == current_phase
    void SwitchMeasuringFromTo(const Phase & current_phase, const Phase & new_phase);

    // TODO: Export total time
    void ExportTimes(void) const;
    void ExportPersons(const std::vector<B18TrafficPerson> & traffic_persons) const;

private:
    Phase current_phase_;
    std::map<Phase, QTime> phases_timers_;
    std::map<Phase, int> phases_total_times_;
};
