import React from 'react';
import ReactWordcloud, {OptionsProp} from 'react-wordcloud';
// import DemographicModal from "../DemographicModal/DemographicModal";
import 'tippy.js/dist/tippy.css';
import 'tippy.js/animations/scale.css';
import {ThemeContext} from "../../contexts/ThemeContext";
// import {isMobile} from "react-device-detect";
// import { Results } from "../../types";
// import {searchTypes} from "../../constants";
import { useTermData } from "../../hooks/useTermData";

const options: OptionsProp  = {
    enableTooltip: true,
    enableOptimizations: true,
    deterministic: true,
    fontFamily: "trebuchet ms",
    fontSizes: [10, 60],
    fontStyle: "normal",
    fontWeight: "normal",
    padding: 1,
    rotations: 0,
    scale: "log",
    spiral: "archimedean",
    transitionDuration: 1000,
}

const WordCloudChart = () => {
    const {theme} = React.useContext(ThemeContext);
    // const [showDemographicModal, setShowDemographicModal] = useState<boolean>(false);
    // const [selectedWord, setSelectedWord] = useState<string | null>(null);

    const { data } = useTermData();

    if (!data) {
        return null;
    }

    const words = data.map((item) => {
        return {
            text: item.x,
            value: item.y,
        }
    })

    const callbacks = {
        // Fix this later type
        getWordColor: (word: { text: string; value: number; }) => {
            const maxFrequency = Math.max(...words.map(w => w.value))
            const frequencyRatio = word.value / maxFrequency

            if (theme === 'dark') {
                const darkGray = [172, 181, 189];
                const brightRed = [255, 0, 0];

                const interpolatedColor = darkGray.map((start, i) => {
                    const end = brightRed[i];
                    if (end !== undefined) {
                        return Math.floor(start + frequencyRatio * (end - start));
                    } else {
                        return start;
                    }
                });

                return `rgb(${interpolatedColor.join(', ')})`;
            } else {
                const redComponent = Math.floor(frequencyRatio * 255);
                return `rgb(${redComponent}, 0, 0)`;
            }
        },
        // Fix this later type
        // onWordClick: (word: { text: string; }) => {
        //     if(props.searchOptions.searchBy === searchTypes[2]) {
        //         setSelectedWord(word.text);
        //         setShowDemographicModal(true);
        //     }
        // }
    };

    // useEffect(() => {
    //     if(!showDemographicModal) {
    //         setSelectedWord(null);
    //     }
    // }, [showDemographicModal])

    return (
        <React.Fragment>
            <ReactWordcloud words={words} options={options} callbacks={callbacks} />
            {/*{props.searchOptions.searchBy === searchTypes[2] &&*/}
            {/*    <DemographicModal*/}
            {/*    show={showDemographicModal}*/}
            {/*    handleClose={() => setShowDemographicModal(false)}*/}
            {/*    selectedWord={selectedWord}*/}
            {/*    isMobile={isMobile}*/}
            {/*/>}*/}
        </React.Fragment>
    )
}

export default WordCloudChart;