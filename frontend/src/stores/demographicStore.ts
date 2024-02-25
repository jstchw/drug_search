import { create } from 'zustand';

interface StoreStates {
  showDemographic: boolean;
  demographicTerm: string;
  demographicType: string;
  groupPageKeys: string[];
}

interface StoreActions {
  setShowDemographic: (showDemographic: boolean) => void;
  setDemographicTerm: (term: string) => void;
  setDemographicType: (searchedType: string) => void;
  setGroupPageKeys: (groupPageKeys: string[]) => void;
}

type Store = StoreStates & StoreActions;

const useDemographicStore = create<Store>((set) => ({
  showDemographic: false,
  setShowDemographic: (showDemographic: boolean) => set({ showDemographic }),
  demographicTerm: '',
  setDemographicTerm: (term: string) => set({ demographicTerm: term }),
  demographicType: '',
  setDemographicType: (searchedType: string) => set({ demographicType: searchedType }),
  groupPageKeys: [],
  setGroupPageKeys: (groupPageKeys: string[]) => set({ groupPageKeys: groupPageKeys }),
}));

export default useDemographicStore;
