import React from 'react'
import {Badge, Col, Modal, OverlayTrigger, Popover, Row} from 'react-bootstrap'
import {useEffect, useState} from "react";
import DrugDescription from "../DrugDescription/DrugDescription";
import DrugAccordion from "../DrugAccordion/DrugAccordion";
import useDrugInfo from "../../hooks/useDrugInfo";
import ApexChart from "../ApexChart/ApexChart";
import { getSideEffectsForDemographics } from "../../services/FDA_Request";
import { searchTypes } from "../OptionModal/OptionModal";

// Function to split the array in chunks,
const chunkArray = (array, size) => {
    const chunked_arr = [];
    let index = 0;
    while (index < array.length) {
        chunked_arr.push(array.slice(index, size + index));
        index += size;
    }
    return chunked_arr;
};

const DemographicModalInfo = ({word}) => {
    const wordRef = React.useRef(null)
    const drugInfo = useDrugInfo([word], searchTypes[0].value)
    const [demographicInfo, setDemographicInfo] = useState([])
    const [isLoading, setIsLoading] = useState(false);

    // Get demographic info only if the word has changed
    useEffect(() => {
        if (word !== wordRef.current) {
            wordRef.current = word;
            setIsLoading(true);
            getSideEffectsForDemographics(word)
                .then((res) => {
                    setDemographicInfo(res);
                })
                .catch((err) => {
                    console.log(err);
                })
                .finally(() => {
                    setIsLoading(false);
                });
        }
    }, [word]);


    const getDemographicDef = (name, count, def) => {
        const popover = (
            <Popover>
                <Popover.Body>
                    <div className={'fs-6 text-center'}>{def}</div>
                    <div className={'fs-6 text-center'}><strong>{count.toLocaleString()}</strong> cases recorded</div>
                </Popover.Body>
            </Popover>
        )

        return (
            <OverlayTrigger trigger={['hover', 'focus']} placement={'top'} overlay={popover}>
                <Badge className={'mx-1 no-pointer-badge'}>{name}</Badge>
            </OverlayTrigger>
        )
    }


    return (
        <div>
            {drugInfo && drugInfo.length > 0 &&
                <React.Fragment>
                    <Row>
                    <Col xs={4}>
                        <DrugDescription
                            drugInfo={drugInfo[0]}
                            style={{border: "none"}}
                        />
                    </Col>
                    <Col xs={4}>
                        <DrugAccordion drugInfo={drugInfo[0]}/>
                    </Col>
                    </Row>
                    <Col className={'text-center'}>
                        {isLoading ?
                            <p>LOADING...</p>
                            :
                            (demographicInfo && Object.keys(demographicInfo).length > 0) ?
                                chunkArray(Object.entries(demographicInfo), 2).map((chunk, index) => (
                                    <Row key={index}>
                                        {chunk.map(([key, value], subIndex) => (
                                            <Col className={'text-center'} key={subIndex}>
                                                <h4>{getDemographicDef(key, value.totalCount, value.def)}</h4>
                                                <ApexChart
                                                    eventDict={value.termCountDict}
                                                    totalCount={value.totalCount}
                                                    type={'searched_group'} />
                                            </Col>
                                        ))}
                                    </Row>
                                ))
                                :
                                <p>No demographic info available</p>
                        }
                    </Col>
                </React.Fragment>
            }
        </div>
    )
}

const DemographicModal = ({ show, handleClose, selectedWord}) => {
    const [word, setWord] = useState(null)

    // Make the word not all caps and check if the word has changed
    useEffect(() => {
        if (selectedWord) {
            const newWord = selectedWord.charAt(0).toUpperCase() + selectedWord.slice(1).toLowerCase();
            setWord(prevWord => newWord !== prevWord ? newWord : prevWord);
        }
    }, [selectedWord]);

    const handleOnHide = () => {
        handleClose();
        setWord(null);
    }


    return (
    <Modal centered show={show} onHide={handleOnHide} size={"xl"}>
      <Modal.Header closeButton>
          {word && <Modal.Title>{word}</Modal.Title>}
      </Modal.Header>
      <Modal.Body>
            {word && <DemographicModalInfo word={word}/>}
      </Modal.Body>
    </Modal>
    )
}

export default DemographicModal