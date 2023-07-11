import React, {useEffect, useState} from 'react';
import ReactWordcloud from 'react-wordcloud';
import DemographicModal from "../DemographicModal/DemographicModal";
import 'tippy.js/dist/tippy.css';
import 'tippy.js/animations/scale.css';
import {searchTypes} from "../OptionModal/OptionModal";

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

    const [showDemographicModal, setShowDemographicModal] = useState(false);
    const [selectedWord, setSelectedWord] = useState(null);

    const data = Object.entries(props.ADEArray)
        .map(([text, value]) => ({text, value}))
        .slice(0, 50)

    const callbacks = {
        getWordColor: word => {
            const maxFrequency = Math.max(...data.map(w => w.value)); // get the maximum frequency
            const frequencyRatio = word.value / maxFrequency; // calculate the ratio of the word's frequency to the maximum frequency
            const redComponent = Math.floor(frequencyRatio * 255); // calculate the red component (0-255)
            return `rgb(${redComponent}, 0, 0)`; // return as an RGB color string
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