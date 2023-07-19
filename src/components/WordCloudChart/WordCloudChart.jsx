import React, {useEffect, useState} from 'react';
import ReactWordcloud from 'react-wordcloud';
import DemographicModal from "../DemographicModal/DemographicModal";
import 'tippy.js/dist/tippy.css';
import 'tippy.js/animations/scale.css';
import {searchTypes} from "../OptionModal/OptionModal";
import {ThemeContext} from "../../contexts/ThemeContext";

const options = {
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
};

const WordCloudChart = (props) => {
    const {theme} = React.useContext(ThemeContext);
    const [showDemographicModal, setShowDemographicModal] = useState(false);
    const [selectedWord, setSelectedWord] = useState(null);

    const data = Object.entries(props.ADEArray)
        .map(([text, value]) => ({text, value}))
        .slice(0, 50)

    const callbacks = {
        // It needs a thorough rework, but it's a start
        getWordColor: word => {
            const maxFrequency = Math.max(...data.map(w => w.value))
            const frequencyRatio = word.value / maxFrequency
            const redComponent = Math.floor(frequencyRatio * 255)
            const greenComponent = theme === 'dark' ? 255 : 0;
            const blueComponent = theme === 'dark' ? Math.floor((1 - frequencyRatio) * 255) : 0;
            return `rgb(${redComponent}, 0, ${blueComponent})`
        },
        onWordClick: (word) => {
            if(props.searchOptions.searchBy === searchTypes[2].value) {
                setSelectedWord(word.text);
                setShowDemographicModal(true);
            }
        }
    };

    useEffect(() => {
        if(!showDemographicModal) {
            setSelectedWord(null);
        }
    }, [showDemographicModal])

    return (
        <React.Fragment>
            <ReactWordcloud words={data} options={options} callbacks={callbacks} />
            {props.searchOptions.searchBy === searchTypes[2].value &&
                <DemographicModal
                show={showDemographicModal}
                handleClose={() => setShowDemographicModal(false)}
                selectedWord={selectedWord}
            />}
        </React.Fragment>
    )
}

export default WordCloudChart;