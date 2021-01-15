from malaria.interventions.malaria_drugs import drug_configs_from_code
from malaria.interventions.malaria_diagnostic import add_diagnostic_survey
from dtk.interventions.triggered_campaign_delay_event import triggered_campaign_delay_event
from dtk.utils.Campaign.CampaignClass import OutbreakIndividualMalaria
from copy import deepcopy, copy
import random
from dtk.utils.Campaign.CampaignClass import *

def add_var_gene_outbreak(cb, start_days: list = None, coverage: float = 1.0,
            repetitions: int = 1, tsteps_btwn_repetitions: int = 60,
            nodeIDs: list = None,
            node_property_restrictions: list = None,
            ind_property_restrictions: list = None,
            irbc_type: list = None,
            msp_type: int = 0,
            minor_epitope_type: list = None
            ):
    """
        Add a var gene specified outbreak
    Args:
        cb: The :py:class:`DTKConfigBuilder <dtk.utils.core.DTKConfigBuilder>`
        object for building, modifying, and writing campaign configuration files.
        start_days: List of integers.
        coverage: Demographic coverage of mda's.
        drug_configs: List of dictionaries of drug configurations to be
            given out, created in add_drug_campaign.
        receiving_drugs_event: (Optional) Broadcast event container with event
            to be broadcast when drugs received.
        repetitions: Number of repetitions for mda. For triggered mda, this is
            for a repeated mda after a trigger.
        tsteps_btwn_repetitions: Timesteps between repeated scheduled mdas or
            between once-triggered repeated mdas.
        nodeIDs: The list of nodes to apply this intervention to (**Node_List**
            parameter). If not provided, set value of NodeSetAll.
        expire_recent_drugs: PropertyValueChanger intervention that updates
            DrugStatus:Recent drug to individual properties.
        node_property_restrictions: List of NodeProperty key:value pairs that nodes
            must have to receive the diagnostic intervention. For example,
            ``[{"NodeProperty1":"PropertyValue1"},
            {"NodeProperty2":"PropertyValue2"}]``. Default is no restrictions.
        disqualifying_properties: List of IndividualProperty key:value pairs that
            cause an intervention to be aborted. For example,
            ``[{"IndividualProperty1":"PropertyValue1"},
            {"IndividualProperty2":"PropertyValue2"}]``.
        target_group: A dictionary targeting an age range and gender of
            individuals for treatment. In the format
            ``{"agemin": x, "agemax": y, "gender": z}``. Default is Everyone.
        trigger_condition_list: List of event triggers upon which mda(s) are
            distributed.
        listening_duration: Duration to listen for the trigger. Default is -1,
            which listens indefinitely.
        triggered_campaign_delay: Delay period between the trigger and the mda.
            Default is 0.
        target_residents_only: When set to true (1), the intervention is only
            distributed to individuals for whom the node is their home node.
            They are not visitors from another node.
        check_eligibility_at_trigger: If triggered event is delayed, you have an
            option to check individual/node's eligibility at the initial trigger
            or when the event is actually distributed after delay.

    Returns:
        None
    """

    if start_days is None:
        start_days = [0]

    if node_property_restrictions is None:
        node_property_restrictions = []
    if ind_property_restrictions is None:
        ind_property_restrictions = []

    if nodeIDs:
        nodeset_config = NodeSetNodeList(Node_List=nodeIDs)
    else:
        nodeset_config = NodeSetAll()


    for start_day in start_days:
        outbreak_event = CampaignEvent(
            Start_Day=start_day,
            Nodeset_Config=nodeset_config,
            Event_Coordinator_Config=StandardInterventionDistributionEventCoordinator(

                Node_Property_Restrictions=node_property_restrictions,
                Property_Restrictions_Within_Node=ind_property_restrictions,
                Demographic_Coverage=coverage,
                Intervention_Config=OutbreakIndividualMalariaVarGenes(
                    IRBC_Type=irbc_type,
                    Ignore_Immunity=True,
                    MSP_Type=msp_type,
                    Minor_Epitope_Type=minor_epitope_type
                    ),

                Number_Repetitions=repetitions,
                Timesteps_Between_Repetitions=tsteps_btwn_repetitions))

        cb.add_event(outbreak_event)


def add_malaria_outbreak(cb, start_days: list = None, coverage: float = 1.0,
            repetitions: int = 1, tsteps_btwn_repetitions: int = 60,
            nodeIDs: list = None,
            node_property_restrictions: list = None,
            ind_property_restrictions: list = None,
            Ignore_Immunity = False,
         ):
    """
        Add a malaria outbreak
    Args:
        cb: The :py:class:`DTKConfigBuilder <dtk.utils.core.DTKConfigBuilder>`
        object for building, modifying, and writing campaign configuration files.
        start_days: List of integers.
        coverage: Demographic coverage of mda's.
        drug_configs: List of dictionaries of drug configurations to be
            given out, created in add_drug_campaign.
        receiving_drugs_event: (Optional) Broadcast event container with event
            to be broadcast when drugs received.
        repetitions: Number of repetitions for mda. For triggered mda, this is
            for a repeated mda after a trigger.
        tsteps_btwn_repetitions: Timesteps between repeated scheduled mdas or
            between once-triggered repeated mdas.
        nodeIDs: The list of nodes to apply this intervention to (**Node_List**
            parameter). If not provided, set value of NodeSetAll.
        expire_recent_drugs: PropertyValueChanger intervention that updates
            DrugStatus:Recent drug to individual properties.
        node_property_restrictions: List of NodeProperty key:value pairs that nodes
            must have to receive the diagnostic intervention. For example,
            ``[{"NodeProperty1":"PropertyValue1"},
            {"NodeProperty2":"PropertyValue2"}]``. Default is no restrictions.
        disqualifying_properties: List of IndividualProperty key:value pairs that
            cause an intervention to be aborted. For example,
            ``[{"IndividualProperty1":"PropertyValue1"},
            {"IndividualProperty2":"PropertyValue2"}]``.
        target_group: A dictionary targeting an age range and gender of
            individuals for treatment. In the format
            ``{"agemin": x, "agemax": y, "gender": z}``. Default is Everyone.
        trigger_condition_list: List of event triggers upon which mda(s) are
            distributed.
        listening_duration: Duration to listen for the trigger. Default is -1,
            which listens indefinitely.
        triggered_campaign_delay: Delay period between the trigger and the mda.
            Default is 0.
        target_residents_only: When set to true (1), the intervention is only
            distributed to individuals for whom the node is their home node.
            They are not visitors from another node.
        check_eligibility_at_trigger: If triggered event is delayed, you have an
            option to check individual/node's eligibility at the initial trigger
            or when the event is actually distributed after delay.

    Returns:
        None
    """

    if start_days is None:
        start_days = [0]

    if node_property_restrictions is None:
        node_property_restrictions = []
    if ind_property_restrictions is None:
        ind_property_restrictions = []

    if nodeIDs:
        nodeset_config = NodeSetNodeList(Node_List=nodeIDs)
    else:
        nodeset_config = NodeSetAll()


    for start_day in start_days:
        outbreak_event = CampaignEvent(
            Start_Day=start_day,
            Nodeset_Config=nodeset_config,
            Event_Coordinator_Config=StandardInterventionDistributionEventCoordinator(

                Node_Property_Restrictions=node_property_restrictions,
                Property_Restrictions_Within_Node=ind_property_restrictions,
                Demographic_Coverage=coverage,
                Intervention_Config=OutbreakIndividualMalaria(

                    Ignore_Immunity=Ignore_Immunity,

                    ),

                Number_Repetitions=repetitions,
                Timesteps_Between_Repetitions=tsteps_btwn_repetitions))

        cb.add_event(outbreak_event)


