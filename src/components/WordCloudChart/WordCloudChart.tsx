// import React from 'react';
// // import DemographicModal from "../DemographicModal/DemographicModal";
// import 'tippy.js/dist/tippy.css';
// import 'tippy.js/animations/scale.css';
// import {ThemeContext} from "../../contexts/ThemeContext";
// // import {isMobile} from "react-device-detect";
// // import { Results } from "../../types";
// // import {searchTypes} from "../../constants";

// const WordCloudChart = () => {
// const {theme} = React.useContext(ThemeContext);
// const [showDemographicModal, setShowDemographicModal] = useState<boolean>(false);
// const [selectedWord, setSelectedWord] = useState<string | null>(null);

// useEffect(() => {
//     if(!showDemographicModal) {
//         setSelectedWord(null);
//     }
// }, [showDemographicModal])

//     return (
//         <React.Fragment>
//             <ReactWordcloud words={words} options={options} callbacks={callbacks} />
//             {/*{props.searchOptions.searchBy === searchTypes[2] &&*/}
//             {/*    <DemographicModal*/}
//             {/*    show={showDemographicModal}*/}
//             {/*    handleClose={() => setShowDemographicModal(false)}*/}
//             {/*    selectedWord={selectedWord}*/}
//             {/*    isMobile={isMobile}*/}
//             {/*/>}*/}
//         </React.Fragment>
//     )
// }

// export default WordCloudChart;
