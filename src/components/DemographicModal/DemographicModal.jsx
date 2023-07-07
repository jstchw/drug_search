import React from 'react'
import { Modal } from 'react-bootstrap'
import {useEffect, useState} from "react";
import DrugDescription from "../DrugDescription/DrugDescription";
import DrugAccordion from "../DrugAccordion/DrugAccordion";
import useDrugInfo from "../../hooks/useDrugInfo";
import {searchTypes} from "../OptionModal/OptionModal";

const DemographicModalInfo = ({word, show}) => {
    const drugInfo = useDrugInfo([word], searchTypes[0].value)

    return (
        <div>
            {drugInfo && drugInfo.length > 0 &&
                <React.Fragment>
                    <DrugDescription drugInfo={drugInfo[0]} showAdditionalSearch={false}/>
                    <DrugAccordion drugInfo={drugInfo[0]}/>
                </React.Fragment>
            }
        </div>
    )
}

const DemographicModal = ({ show, handleClose, selectedWord}) => {
    const [word, setWord] = useState(null)

    // Make the word not all caps
    useEffect(() => {
        if (selectedWord) {
            const newWord = selectedWord.charAt(0).toUpperCase() + selectedWord.slice(1).toLowerCase()
            if(newWord !== word) {
                setWord(newWord)
            }
        }
    }, [selectedWord])

    const handleModalClose = () => {
        handleClose()
        setWord(null)
    }



    return (
    <Modal centered show={show} onHide={handleModalClose}>
      <Modal.Header closeButton>
          {word && <Modal.Title>{word}</Modal.Title>}
      </Modal.Header>
      <Modal.Body>
            {word && <DemographicModalInfo word={word} show={show}/>}
      </Modal.Body>
    </Modal>
    )
}

export default DemographicModal